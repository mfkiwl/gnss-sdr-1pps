/*!
 * \file 
 * \brief  Interface of a GNU Radio block that sends a STOP message to the
 * control queue right after a specific number of samples have passed through it.
 * \author Javier Arribas, 2018. jarribas(at)cttc.es
 * \author Carlos Aviles, 2010. carlos.avilesr(at)googlemail.com
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2019  (see AUTHORS file for a list of contributors)
 *
 * GNSS-SDR is a software defined Global Navigation
 *          Satellite Systems receiver
 *
 * This file is part of GNSS-SDR.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * -------------------------------------------------------------------------
 */


#ifndef GNSS_SDR_GNSS_SDR_SPOOF_H
#define GNSS_SDR_GNSS_SDR_SPOOF_H

#include "concurrent_queue.h"
#include <boost/shared_ptr.hpp>
#include <gnuradio/sync_block.h>  // for sync_block
#include <gnuradio/types.h>       // for gr_vector_const_void_star
#include <pmt/pmt.h>
#include <cstddef>  // for size_t
#include <cstdint>
#include <memory>

#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/fft_shift.h>

#define CHUNK_SIZE (2048*8*2) // ~1023 MS/s/32768=30~Hz/bin
#define KEEP_SIZE   (600)     // 30 Hz/bin * 600 = ~+/-20 kHz

class Gnss_Spoofing_Protect;

boost::shared_ptr<Gnss_Spoofing_Protect> gnss_sdr_make_spoof(
    size_t sizeof_stream_item,
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue);

boost::shared_ptr<Gnss_Spoofing_Protect> gnss_sdr_make_spoof(
    size_t sizeof_stream_item,
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue);

/*!
 * \brief Implementation of a GNU Radio block that sends a STOP message to the
 * control queue right after a specific number of samples have passed through it.
 */
class Gnss_Spoofing_Protect : public gr::sync_block
{
public:
    int work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items);

private:
    friend boost::shared_ptr<Gnss_Spoofing_Protect> gnss_sdr_make_spoof(
        size_t sizeof_stream_item,
        std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue);

    Gnss_Spoofing_Protect(size_t sizeof_stream_item,
        std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue);

    gr::fft::fft_complex* plan = new gr::fft::fft_complex(CHUNK_SIZE, true);
    gr_complex spoofing_average_mul[KEEP_SIZE*2]; // 2* since sta followed by sto -> must fftshift to put 0 at center
    gr_complex spoofing_average_div[KEEP_SIZE*2]; // 2* since sta followed by sto -> must fftshift to put 0 at center
    gr_complex bufout0_sta[KEEP_SIZE];
    gr_complex bufout0_sto[KEEP_SIZE];
    gr_complex bufout_sta[KEEP_SIZE];
    gr_complex bufout_sto[KEEP_SIZE];
//    gr_complex bufout_tmpsta[KEEP_SIZE];
//    gr_complex bufout_tmpsto[KEEP_SIZE];
    gr_complex processed_output[CHUNK_SIZE];
    gr_complex weight_;
    float stdargres_;
    int avg_index_;
    int num_file_;
    uint64_t d_ncopied_items;
    int first_time_;
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> d_queue;
};

#endif  // GNSS_SDR_GNSS_SDR_SPOOF_H
