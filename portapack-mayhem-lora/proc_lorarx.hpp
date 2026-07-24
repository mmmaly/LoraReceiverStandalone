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

#ifndef __PROC_LORARX_H__
#define __PROC_LORARX_H__

#include "baseband_processor.hpp"
#include "baseband_thread.hpp"
#include "rssi_thread.hpp"
#include "message.hpp"
#include "dsp_decimate.hpp"

// The portable, no-malloc LoRa demodulator core (verified on real captures).
#define LORA_LITE_OS 4
#define LORA_LITE_MAX_SF 8
#define LORA_LITE_MAX_PAYLOAD 64
#include "lora_lite.h"

// Baseband: SR = bw * os * 8. For bw = 62.5 kHz -> 2.0 MHz in, /8 -> 250 kHz
// (os = 4). Two-stage FIR decimation (same filters as the weather proc).
class LoRaRxProcessor : public BasebandProcessor {
   public:
    void execute(const buffer_c8_t& buffer) override;
    void on_message(const Message* const message) override;

   private:
    static void packet_callback(const lora_lite::Packet& p, void* user);
    void configure(const LoRaRxConfigureMessage& message);

    bool configured{false};
    size_t baseband_fs = 0;

    // /4 then /2 = /8 total, identical to proc_weather
    std::array<complex16_t, 512> dst{};
    const buffer_c16_t dst_buffer{dst.data(), dst.size()};
    dsp::decimate::FIRC8xR16x24FS4Decim4 decim_0{};
    dsp::decimate::FIRC16xR16x16Decim2 decim_1{};

    // Scratch for the float conversion handed to the core.
    std::array<lora_lite::Cx, 256> cx_buf{};

    // NB: the ~63 KB lora_lite::Demod is NOT a member. This processor is
    // heap-allocated (make_unique) and the M4 heap can't hold 63 KB, so the
    // demod lives in static storage (.bss) instead. See proc_lorarx.cpp.

    BasebandThread baseband_thread{baseband_fs, this, baseband::Direction::Receive};
    RSSIThread rssi_thread{};
};

#endif /*__PROC_LORARX_H__*/
