/* -*- c++ -*- */
/*
 * Copyright 2020 gr-sgd author.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_SGD_SGD_IMPL_H
#define INCLUDED_SGD_SGD_IMPL_H
#if GNURADIO_USES_STD_POINTERS
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif
#include <gnuradio/sync_block.h>  // for sync_block
#include <gnuradio/types.h>       // for gr_vector_const_void_star
#include <pmt/pmt.h>

#include <volk/volk.h>

class sgd_impl;

#if GNURADIO_USES_STD_POINTERS
std::shared_ptr<sgd_impl> gnss_sdr_make_sgd(
    int delay_max, float seuil, float alpha,
        bool mean, int mean_length, int iter_count);
#else
boost::shared_ptr<sgd_impl> gnss_sdr_make_sgd(
    int delay_max, float seuil, float alpha,
        bool mean, int mean_length, int iter_count);
#endif

    class sgd_impl : public gr::sync_block
    {
     private:
#if GNURADIO_USES_STD_POINTERS
         friend std::shared_ptr<sgd_impl> gnss_sdr_make_sgd(
             int delay_max, float seuil, float alpha,
                  bool mean, int mean_length, int iter_count);
#else
         friend boost::shared_ptr<sgd_impl> gnss_sdr_make_sgd(
             int delay_max, float seuil, float alpha,
                  bool mean, int mean_length, int iter_count);
#endif

      int _w1_size;
      float _seuil;
      float _alpha;
      int _delay_max;
      uint32_t alignment;

      gr_complex *tmp2;
      gr_complex e;
      float val_max;
      uint32_t *val_max_index;
      gr_complex *XXxw1;
      gr_complex *w1;
      float *w1_mag;
      gr_complex *xxConj;
      int _residual;
      FILE *_w1_out;
      int _iter;
      std::complex<double> **_w1_array;
      std::complex<double> *w1_accum;
      gr_complex *w1_res;
      int _array_index;
      int _nb_accum;
      bool _mean;
      int _mean_length;
      int _iter_count;

     public:
      sgd_impl();
      sgd_impl(int w1_size, float seuil, float alpha,
                bool mean, int mean_length, int iter_count);
      ~sgd_impl();
      int fixed_rate_ninput_to_noutput(int ninput);

      // Where all the action really happens
      int work(
              int noutput_items,
              gr_vector_const_void_star &input_items,
              gr_vector_void_star &output_items
      );
    };

#endif /* INCLUDED_SGD_SGD_IMPL_H */

