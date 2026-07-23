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

#include "ui.hpp"
#include "ui_lorarx.hpp"
#include "ui_navigation.hpp"
#include "external_app.hpp"

namespace ui::external_app::lorarx {
void initialize_app(ui::NavigationView& nav) {
    nav.push<LoRaRxView>();
}
}  // namespace ui::external_app::lorarx

extern "C" {

__attribute__((section(".external_app.app_lorarx.application_information"), used)) application_information_t _application_information_lorarx = {
    /*.memory_location = */ (uint8_t*)0x00000000,
    /*.externalAppEntry = */ ui::external_app::lorarx::initialize_app,
    /*.header_version = */ CURRENT_HEADER_VERSION,
    /*.app_version = */ VERSION_MD5,

    /*.app_name = */ "LoRa RX",
    /*.bitmap_data = */ {
        0x00,
        0x00,
        0x7C,
        0x3E,
        0x44,
        0x22,
        0x44,
        0x22,
        0x44,
        0x22,
        0x7C,
        0x3E,
        0x10,
        0x08,
        0x10,
        0x08,
        0x10,
        0x08,
        0xFE,
        0x7F,
        0x82,
        0x41,
        0xBA,
        0x5D,
        0xAA,
        0x55,
        0xBA,
        0x5D,
        0x82,
        0x41,
        0xFE,
        0x7F,
    },
    /*.icon_color = */ ui::Color::green().v,
    /*.menu_location = */ app_location_t::RX,
    /*.desired_menu_position = */ -1,

    /*.m4_app_tag = portapack::spi_flash::image_tag_lorarx */ {'P', 'L', 'O', 'R'},
    /*.m4_app_offset = */ 0x00000000,  // will be filled at compile time
};
}
