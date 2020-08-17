/*!
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

#include "jamming_detection.h"
#include <glog/logging.h>           // for LOG
#include <gnuradio/io_signature.h>  // for io_signature
#include <algorithm>                // for min
#include <cstring>                  // for memcpy
#include <unistd.h>                 // for usleep
#include <utility>
#include <volk/volk.h>

#pragma message("Jamming detection compile")

#define MEMORY_LEN 5 // remember jamming even std rises suddendly

Gnss_Jamming_Protect::Gnss_Jamming_Protect(float threshold, int averages) : gr::sync_block("jamming_detection",
                               gr::io_signature::make(1, 20, sizeof(gr_complex)),
                               gr::io_signature::make(1, 1, sizeof(gr_complex))),
                           d_threshold(threshold),
                           d_averages(averages)
{
    printf("Gnss_Jamming_Protect\n");
/*
Ron Economos (April 5, 2020 10:58 AM)
To: discuss-gnuradio@gnu.org
I would use set_output_multiple() instead. See my previous e-mail for an
example.
https://lists.gnu.org/archive/html/discuss-gnuradio/2019-08/msg00188.html
*/
    set_output_multiple(CHUNK_SIZE); // only trigger processing if that amount of samples was accumulated
    avg_index_=0;
    weight_={0.,0.};
    jamming_memory_=0;
}

#if GNURADIO_USES_STD_POINTERS
std::shared_ptr<Gnss_Jamming_Protect> gnss_sdr_make_jamm(float threshold, int averages)
{
    std::shared_ptr<Gnss_Jamming_Protect> jamming_detect_(new Gnss_Jamming_Protect(threshold, averages));
    return jamming_detect_;
}
#else
boost::shared_ptr<Gnss_Jamming_Protect> gnss_sdr_make_jamm(float threshold, int averages)
{
    boost::shared_ptr<Gnss_Jamming_Protect> jamming_detection(new Gnss_Jamming_Protect(threshold, averages));
    printf("Jamming detection: variable created\n");
    return jamming_detection;
}
#endif

int Gnss_Jamming_Protect::work(int noutput_items,
    gr_vector_const_void_star &input_items,
    gr_vector_void_star &output_items)
{   long unsigned int ch = 0;
    unsigned int alignment = volk_get_alignment();
    gr_complex *bufin; // ,*bufin;
    gr_complex integral;
    const gr_complex* in;
    gr_complex* carre=(gr_complex*)volk_malloc(sizeof(gr_complex)*CHUNK_SIZE, alignment);
// see https://github.com/gnss-sdr/gnss-sdr/blob/master/src/algorithms/acquisition/gnuradio_blocks/pcps_acquisition_fine_doppler_cc.h for declaration of gr::fft
    bufin=plan->get_inbuf();
    // ibufin=iplan->get_inbuf();
    for (ch = 0; ch < input_items.size(); ch++)
        { // identity: output the same as 1st channel input
          in= (const gr_complex*)input_items[ch]; // all channels
          memcpy(bufin, in, CHUNK_SIZE * sizeof(gr_complex));
          plan->execute();
          if (ch==0)
             {memcpy(bufout0,plan->get_outbuf(),CHUNK_SIZE * sizeof(gr_complex)); // save FFT(CH0)
             }
          if (ch==1)
             {
              // volk_32fc_x2_multiply_conjugate_32fc(bufout,plan->get_outbuf(),bufout0,CHUNK_SIZE);
              volk_32fc_x2_divide_32fc(bufout,plan->get_outbuf(),bufout0,CHUNK_SIZE); // CH1/CH0
              //memcpy(ibufin, bufout, CHUNK_SIZE * sizeof(gr_complex));
              //iplan->execute(); // result in iplan->get_outbuf()
              //weight_+=iplan->get_outbuf()[0]; 
              integral={0.,0.};
              for (int i=0;i<CHUNK_SIZE;i++) integral+=bufout[i];
              integral/=CHUNK_SIZE;
              weight_avg_+=integral;
              avg_index_++;
             }
        }
    if (avg_index_==d_averages)    // restart averaging
      {weight_=weight_avg_/(float)avg_index_;
       weight_avg_={0.,0.};
       avg_index_=0;
       printf("xcorr: %f+i*%f -> %f ",weight_.real(),weight_.imag(),norm(weight_));
      }
      if ((norm(weight_)<d_threshold)&&(jamming_memory_==0)) // no jamming
       {
        memcpy(output_items[0], input_items[0], noutput_items * input_signature()->sizeof_stream_item(ch));
        printf("\n");
       }
    else
      {if (norm(weight_)>=d_threshold)
          {printf(" /!\\\n");
           jamming_memory_=MEMORY_LEN; // reinit memory
          }
       else
          {printf(" \\!/\n");
           jamming_memory_--;  // |weight_| < NORM_THRESHOLD so jamming_memory_>0 and --
          }
       volk_32fc_s32fc_multiply_32fc(carre,(const gr_complex*)input_items[0],-weight_,CHUNK_SIZE); 
       volk_32fc_x2_add_32fc((gr_complex*)output_items[0], (const gr_complex*)input_items[1], carre, CHUNK_SIZE);
      }
//    delete plan;
    volk_free(carre);
    return noutput_items;
}
