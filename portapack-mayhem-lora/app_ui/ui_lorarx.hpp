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

#ifndef __UI_LORARX_H__
#define __UI_LORARX_H__

#include "ui.hpp"
#include "ui_navigation.hpp"
#include "ui_receiver.hpp"
#include "ui_freq_field.hpp"
#include "app_settings.hpp"
#include "radio_state.hpp"
#include "log_file.hpp"

using namespace ui;

namespace ui::external_app::lorarx {

class LoRaRxView : public View {
   public:
    LoRaRxView(NavigationView& nav);
    ~LoRaRxView();

    void focus() override;
    std::string title() const override { return "LoRa RX"; };

   private:
    void configure_baseband();
    void on_packet(const LoRaPacketData& pkt);

    NavigationView& nav_;
    RxRadioState radio_state_{};
    app_settings::SettingsManager settings_{"rx_lora", app_settings::Mode::RX};

    uint32_t packet_count{0};

    RFAmpField field_rf_amp{{13 * 8, 0 * 16}};
    LNAGainField field_lna{{15 * 8, 0 * 16}};
    VGAGainField field_vga{{18 * 8, 0 * 16}};
    RSSI rssi{{21 * 8, 0, 6 * 8, 4}};
    Channel channel{{21 * 8, 5, 6 * 8, 4}};

    RxFrequencyField field_frequency{{0 * 8, 0 * 16}, nav_};

    // SF 7..12
    OptionsField options_sf{
        {0 * 8, 1 * 16 + 8},
        6,
        {{"SF7", 7}, {"SF8", 8}}};

    // Bandwidth (kHz) — decimation supports 62.5k and 125k
    OptionsField options_bw{
        {7 * 8, 1 * 16 + 8},
        8,
        {{"BW62.5k", 62500}, {"BW125k", 125000}}};

    // Sync word (0x12 public, 0x34 private, 0 = any)
    OptionsField options_sync{
        {17 * 8, 1 * 16 + 8},
        9,
        {{"sync 12", 0x12}, {"sync 34", 0x34}, {"sync any", 0x00}}};

    Console console{{0, 3 * 16, screen_width, screen_height - 3 * 16}};

    MessageHandlerRegistration message_handler_packet{
        Message::ID::LoRaRxPacket,
        [this](Message* const p) {
            const auto message = static_cast<const LoRaRxPacketMessage*>(p);
            this->on_packet(*message->packet);
        }};
};

}  // namespace ui::external_app::lorarx

#endif /*__UI_LORARX_H__*/
