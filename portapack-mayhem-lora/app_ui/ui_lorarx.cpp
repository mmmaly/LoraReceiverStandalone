/*
 * Copyright (C) 2026 (LoRa RX app)
 *
 * This file is part of PortaPack.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include "ui_lorarx.hpp"
#include "baseband_api.hpp"
#include "string_format.hpp"
#include "portapack.hpp"

using namespace portapack;

namespace ui::external_app::lorarx {

LoRaRxView::LoRaRxView(NavigationView& nav)
    : nav_{nav} {
    baseband::run_prepared_image(portapack::memory::map::m4_code.base());

    add_children({&field_rf_amp,
                  &field_lna,
                  &field_vga,
                  &rssi,
                  &channel,
                  &field_frequency,
                  &options_sf,
                  &options_bw,
                  &options_sync,
                  &console});

    field_frequency.set_value(869618000);  // Slovak MeshCore SF8 channel
    options_sf.set_selected_index(1);      // SF8
    options_bw.set_selected_index(0);      // 62.5 kHz
    options_sync.set_selected_index(0);    // 0x12

    options_sf.on_change = [this](size_t, auto) { configure_baseband(); };
    options_bw.on_change = [this](size_t, auto) { configure_baseband(); };
    options_sync.on_change = [this](size_t, auto) { configure_baseband(); };

    receiver_model.enable();
    configure_baseband();
    console.writeln("LoRa RX ready. Waiting for packets...");
}

void LoRaRxView::configure_baseband() {
    const uint8_t sf = (uint8_t)options_sf.selected_index_value();
    const uint32_t bw = (uint32_t)options_bw.selected_index_value();
    const uint16_t sync = (uint16_t)options_sync.selected_index_value();

    // Baseband sample rate = bw * os(4) * decim(8). Keep the radio in step.
    const uint32_t fs = bw * 4 * 8;
    receiver_model.set_sampling_rate(fs);
    receiver_model.set_baseband_bandwidth(fs / 2);

    // cr is taken from the explicit header; pass 1 as a default, CRC on.
    baseband::set_lorarx(sf, 1, true, sync, bw);
}

void LoRaRxView::on_packet(const LoRaPacketData& pkt) {
    packet_count++;

    // Console colour: ESC followed by a palette index (see afsk_rx). 11=green,
    // 9=red for CRC ok / bad.
    std::string line = "\x1B";
    line += (char)(pkt.crc_ok ? 11 : 9);
    line += to_string_dec_uint(packet_count) + ": ";
    for (uint8_t i = 0; i < pkt.len; i++)
        line += to_string_hex(pkt.data[i], 2);
    console.writeln(line);

    std::string info = "  SNR " + to_string_dec_int((int)pkt.snr) +
                       " sync " + to_string_hex(pkt.sync, 2) +
                       (pkt.crc_ok ? " CRC ok" : " CRC BAD");
    console.writeln(info);
}

void LoRaRxView::focus() {
    field_frequency.focus();
}

LoRaRxView::~LoRaRxView() {
    receiver_model.disable();
    baseband::shutdown();
}

}  // namespace ui::external_app::lorarx
