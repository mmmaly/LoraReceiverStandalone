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

#include "proc_lorarx.hpp"
#include "portapack_shared_memory.hpp"
#include "event_m4.hpp"

#include <cstddef>

// The demodulator (~63 KB of buffers) lives in static storage, not on the M4
// heap. Placing it in .bss keeps the heap-allocated processor object tiny so
// its construction can't fail for lack of heap.
static lora_lite::Demod g_demod;

void LoRaRxProcessor::execute(const buffer_c8_t& buffer) {
    if (!configured) return;

    // 2048 complex8 in -> /4 -> 512 -> /2 -> 256 complex16 out (os=4 @ bw)
    const auto decim_0_out = decim_0.execute(buffer, dst_buffer);
    const auto decim_1_out = decim_1.execute(decim_0_out, dst_buffer);
    feed_channel_stats(decim_1_out);

    const size_t n = decim_1_out.count;
    for (size_t i = 0; i < n && i < cx_buf.size(); i++) {
        // complex16 -> normalized float, matching the host core's scaling
        cx_buf[i].r = (float)decim_1_out.p[i].real() * (1.0f / 32768.0f);
        cx_buf[i].i = (float)decim_1_out.p[i].imag() * (1.0f / 32768.0f);
    }
    g_demod.feed(cx_buf.data(), (int)(n < cx_buf.size() ? n : cx_buf.size()));
}

void LoRaRxProcessor::packet_callback(const lora_lite::Packet& p, void* /*user*/) {
    auto& pd = *reinterpret_cast<LoRaPacketData*>(shared_memory.bb_data.data);
    pd.len = p.len;
    pd.crc_ok = p.crc_ok ? 1 : 0;
    pd.has_crc = p.has_crc ? 1 : 0;
    pd.sync = p.sync;
    pd.snr = p.snr;
    pd.cfo = p.cfo;
    for (uint8_t i = 0; i < p.len && i < sizeof(pd.data); i++) pd.data[i] = p.data[i];

    LoRaRxPacketMessage message{&pd};
    shared_memory.application_queue.push(message);
}

void LoRaRxProcessor::on_message(const Message* const message) {
    switch (message->id) {
        case Message::ID::LoRaRxConfigure:
            configure(*reinterpret_cast<const LoRaRxConfigureMessage*>(message));
            break;
        default:
            break;
    }
}

void LoRaRxProcessor::configure(const LoRaRxConfigureMessage& message) {
    baseband_fs = message.bandwidth * LORA_LITE_OS * 8;  // /8 decimation chain
    baseband_thread.set_sampling_rate(baseband_fs);

    decim_0.configure(taps_200k_wfm_decim_0.taps);
    decim_1.configure(taps_200k_wfm_decim_1.taps);

    g_demod.init(message.sf, message.cr, message.has_crc != 0, message.sync_word,
               &LoRaRxProcessor::packet_callback, this,
               message.bandwidth, 869618000);

    configured = true;
}

int main() {
    EventDispatcher event_dispatcher{std::make_unique<LoRaRxProcessor>()};
    event_dispatcher.run();
    return 0;
}
