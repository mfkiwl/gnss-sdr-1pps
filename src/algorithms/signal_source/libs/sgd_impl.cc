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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "sgd.h"

//#define DEBUG_COEFF

#if GNURADIO_USES_STD_POINTERS
std::shared_ptr<sgd_impl> gnss_sdr_make_sgd(
    int delay_max, float seuil, float alpha)
{
	std::shared_ptr<sgd_impl> sgd_inst(new sgd_impl(delay_max, seuil, alpha));
	return sgd_inst;
}
#else
boost::shared_ptr<sgd_impl> gnss_sdr_make_sgd(
    int delay_max, float seuil, float alpha)
{
	boost::shared_ptr<sgd_impl> sgd_inst(new sgd_impl(delay_max, seuil, alpha));
	return sgd_inst;
}
#endif



    /*sgd::sptr
    sgd::make(int delay_max, float seuil, float alpha)
    {
      return gnuradio::get_initial_sptr
        (new sgd_impl(delay_max, seuil, alpha));
    }*/


    /*
     * The private constructor
     */
    sgd_impl::sgd_impl(int delay_max, float seuil, float alpha)
      : gr::sync_block("sgd",
              gr::io_signature::make(2, 2, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex))),
              _w1_size(delay_max * 2 + 1), _seuil(seuil), _alpha(alpha),
              _delay_max(delay_max), _residual(0), _w1_out(NULL), _iter(1)
    {
        alignment = volk_get_alignment();
        tmp2   = (gr_complex *)volk_malloc(sizeof(gr_complex) * _w1_size, alignment);
        val_max_index = (uint32_t *)volk_malloc(sizeof(uint32_t), alignment);
        XXxw1  = (gr_complex *)volk_malloc(sizeof(gr_complex), alignment);
        w1     = (gr_complex *)volk_malloc(sizeof(gr_complex) * _w1_size, alignment);
        w1_mag = (float *)volk_malloc(sizeof(float)* _w1_size, alignment);
        xxConj = (gr_complex *)volk_malloc(sizeof(gr_complex) * _w1_size, alignment);

        /* cleanup w1 */
        for (int i=0; i < _w1_size; i++)
            w1[i] = 0;

//#ifdef DEBUG_COEFF
	_w1_out = fopen("w1_dump.bin", "w+");
//#endif
	}

    /*
     * Our virtual destructor.
     */
    sgd_impl::~sgd_impl()
    {
        for (int index = 0; index < _w1_size; index++) {
			float tr = w1[index].real();
			float ti = w1[index].imag();
			fwrite(&tr, sizeof(float), 1, _w1_out);
			fwrite(&ti, sizeof(float), 1, _w1_out);
		}
//#ifdef DEBUG_COEFF
		fclose(_w1_out);
//#endif
        volk_free(tmp2);
        volk_free(val_max_index);
        volk_free(XXxw1);
        volk_free(w1);
        volk_free(w1_mag);
        volk_free(xxConj);
    }

    int
    sgd_impl::fixed_rate_ninput_to_noutput(int ninput)  
    {
      return ninput - (_w1_size);
    }

    int
    sgd_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *x = (const gr_complex *) input_items[0];
      const gr_complex *y = (const gr_complex *) input_items[1];
      gr_complex *res = (gr_complex *) output_items[0];

      int size = noutput_items;

	  if (size < _w1_size) {
		  _residual = size;
		  set_history(_residual);
		  return 0;
	  }

      for (int k=0; k < size-_w1_size; k++) {

        volk_32fc_x2_dot_prod_32fc(XXxw1, x+k, w1, _w1_size);
        //e = (y[k + _delay_max] - *XXxw1) * (float)_alpha;//(1/sqrtf(_iter));
        e = (y[k + _delay_max] - *XXxw1) * (float)(1/sqrtf((float)_iter));
		if (_iter < 10000)
			_iter++;
		else {
			_iter = 1;
			for (int index = 0; index < _w1_size; index ++)
				std::cout << index - (_delay_max) << " " << w1[index] << std::endl;
			//std::cout << " " << std::endl;
		}

        volk_32fc_conjugate_32fc(xxConj, x+k, _w1_size);
        volk_32fc_s32fc_multiply_32fc(tmp2, xxConj, e, _w1_size);

        volk_32fc_x2_add_32fc(w1, w1, tmp2, _w1_size);

// remove if _w1_size == 1 (case of GPS un-jamming)
        /* abs */
        volk_32fc_magnitude_32f(w1_mag, w1, _w1_size);
        /* search max */
        volk_32f_index_max_32u(val_max_index, w1_mag, _w1_size);
        /* store max */
        val_max = w1_mag[*val_max_index];

        val_max *= _seuil;
        /* all entry < val_max => 0 */
        for (int index = 0; index < _w1_size; index++) {
            if (w1_mag[index] < val_max)
                w1[index] = 0;
#ifdef DEBUG_COEFF
			float toto = w1[index].real();
			fwrite(&toto, sizeof(float), 1, _w1_out);
			toto = w1[index].imag();
			fwrite(&toto, sizeof(float), 1, _w1_out);
#endif
		}
      }
// remove until here if _w1_size == 1 (case of GPS un-jamming)

      /* suppress part */
      int maxsize = size;
      int coeffsize = _w1_size;
      int sum_size = maxsize-coeffsize;
      gr_complex tmp;
      for (int i=0; i < sum_size; i++) {
          tmp = 0;
          volk_32fc_x2_dot_prod_32fc(&tmp, x+i, w1, coeffsize);
          res[i] = y[_delay_max + i] - tmp;
      }

      _residual = _w1_size - 1;
	  set_history(_residual);

      // Tell runtime system how many output items we produced.
      return sum_size;
    }


