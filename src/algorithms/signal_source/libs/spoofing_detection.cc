/*!
 * \file gnss_sdr_valve.cc
 * \brief Implementation of a GNU Radio block that sends a STOP message to the
 * control queue right after a specific number of samples have passed through it.
 * \author Javier Arribas, 2018. jarribas(at)cttc.es
 * \author Carlos Aviles, 2010. carlos.avilesr(at)googlemail.com
 *
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

#include "spoofing_detection.h"
#include "command_event.h"
#include <glog/logging.h>           // for LOG
#include <gnuradio/io_signature.h>  // for io_signature
#include <algorithm>                // for min
#include <cstring>                  // for memcpy
#include <unistd.h>                 // for usleep
#include <utility>
#include <volk/volk.h>
//#define CHUNK_SIZE (2048*8*2) // ~1023 MS/s/16384=30~Hz/bin
#define KEEPSIZE 600 // 30 Hz/bin * 600 = ~+/-20 kHz
#define STD_THRESHOLD 0.5 // rad
#define MAXSAT  20  // too many satellites will start detecting genuine constellation ?
#define MAXKEEP 7   // analyze MAXSAT but only keep MAXKEEP

#pragma message("spoofing detection")

#define Navg 16 // FFT averages
#undef avgcpl   // averaging the complex weight rather than its magnitude and phase

Gnss_Spoofing_Protect::Gnss_Spoofing_Protect(size_t sizeof_stream_item,
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue) : gr::sync_block("spoofing_detection",
                               gr::io_signature::make(1, 20, sizeof_stream_item),
                               gr::io_signature::make(1, 1, sizeof_stream_item)),
                           d_ncopied_items(0),
                           d_queue(std::move(queue))
{
    printf("Spoofing Detect\n");
/*
Ron Economos (April 5, 2020 10:58 AM)
To: discuss-gnuradio@gnu.org
I would use set_output_multiple() instead. See my previous e-mail for an 
example.
https://lists.gnu.org/archive/html/discuss-gnuradio/2019-08/msg00188.html
*/
    set_output_multiple(CHUNK_SIZE); // only trigger processing if that amount of samples was accumulated
    first_time_=0;
    memset(spoofing_average_mul,0,sizeof(gr_complex)*KEEP_SIZE*2);
    memset(spoofing_average_div,0,sizeof(gr_complex)*KEEP_SIZE*2);
    avg_index_=0;
    num_file_=0;
    weight_=1.;
    stdargres_=10.;
}


boost::shared_ptr<Gnss_Spoofing_Protect> gnss_sdr_make_spoof(size_t sizeof_stream_item, std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue)
{
//    unsigned int alignment = volk_get_alignment();
    boost::shared_ptr<Gnss_Spoofing_Protect> spoofing_detect_(new Gnss_Spoofing_Protect(sizeof_stream_item, std::move(queue)));
    return spoofing_detect_;
}

int Gnss_Spoofing_Protect::work(int noutput_items,
    gr_vector_const_void_star &input_items,
    gr_vector_void_star &output_items)
{   long unsigned int ch = 0;
    int i,c,cnt;
    uint16_t maxpos;
    unsigned int alignment = volk_get_alignment();
    float maxval,maxvallim,meanarg,meanabs,stdarg[MAXSAT],stdabs[MAXSAT],weightabs=0.,weightarg=0.;
    gr_complex* bufin;
#ifdef avgcpl
    gr_complex stddiv[MAXSAT];
#endif
    gr_complex meandiv;
/*
    FILE *fo;
    char filenam[256];
*/
    int count;
    const gr_complex* in;
    gr_complex* carre=(gr_complex*)volk_malloc(sizeof(gr_complex)*CHUNK_SIZE, alignment);
// see https://github.com/gnss-sdr/gnss-sdr/blob/master/src/algorithms/acquisition/gnuradio_blocks/pcps_acquisition_fine_doppler_cc.h for declaration of gr::fft
//    gr::fft::fft_complex* plan = new gr::fft::fft_complex(CHUNK_SIZE, true);
    bufin=plan->get_inbuf();
//    printf("block: %d items, %ld out, %ld in\n",noutput_items,output_items.size(),input_items.size());
//            block: 32768 items, 1 out, 2 in
    if (input_items.size()!=2) first_time_=1; // don't save if other than 2 input channels
    for (ch = 0; ch < input_items.size(); ch++)
        { // identity: output the same as 1st channel input
          in= (const gr_complex*)input_items[ch]; // all channels
          volk_32fc_x2_multiply_32fc(carre, in, in, CHUNK_SIZE);
          memcpy(bufin, carre, CHUNK_SIZE * sizeof(gr_complex));
          plan->execute();
          if (ch==0) 
             {memcpy(bufout0_sta,plan->get_outbuf(),KEEP_SIZE * sizeof(gr_complex)); // save FFT(CH0)
     //         volk_32fc_s32fc_multiply_32fc(bufout0_sta,bufout0_sta,{1000.0,0.0},KEEP_SIZE);
              memcpy(bufout0_sto,&plan->get_outbuf()[CHUNK_SIZE-KEEP_SIZE-1],KEEP_SIZE * sizeof(gr_complex)); // save FFT(CH0)
     //         volk_32fc_s32fc_multiply_32fc(bufout0_sto,bufout0_sto,{1000.0,0.0},KEEP_SIZE);
             }
          if (ch==1) 
             {// volk_32fc_s32fc_multiply_32fc(bufout_tmpsta,plan->get_outbuf(),{1000.0,0.0},KEEP_SIZE);
              // volk_32fc_s32fc_multiply_32fc(bufout_tmpsto,&plan->get_outbuf()[CHUNK_SIZE-KEEP_SIZE-1],{1000.0,0.0},KEEP_SIZE);
              volk_32fc_x2_multiply_conjugate_32fc(bufout_sta,plan->get_outbuf(),bufout0_sta,KEEP_SIZE); 
              volk_32fc_x2_multiply_conjugate_32fc(bufout_sto,&plan->get_outbuf()[CHUNK_SIZE-KEEP_SIZE-1],bufout0_sto,KEEP_SIZE); 
              volk_32fc_x2_add_32fc( spoofing_average_mul,            spoofing_average_mul,           bufout_sta,KEEP_SIZE); 
              volk_32fc_x2_add_32fc(&spoofing_average_mul[KEEP_SIZE],&spoofing_average_mul[KEEP_SIZE],bufout_sto,KEEP_SIZE); 

              volk_32fc_x2_divide_32fc(bufout_sta,plan->get_outbuf(),bufout0_sta,KEEP_SIZE);  // CH1/CH0
              volk_32fc_x2_divide_32fc(bufout_sto,&plan->get_outbuf()[CHUNK_SIZE-KEEP_SIZE-1],bufout0_sto,KEEP_SIZE); 
              volk_32fc_x2_add_32fc( spoofing_average_div,            spoofing_average_div,           bufout_sta,KEEP_SIZE); 
              volk_32fc_x2_add_32fc(&spoofing_average_div[KEEP_SIZE],&spoofing_average_div[KEEP_SIZE],bufout_sto,KEEP_SIZE); 
              avg_index_++;
             }
/*
	  if ((first_time_==0)&&(avg_index_==Navg))
	     {//if (ch==input_items.size()-1) first_time_=1;                   // save both channels in file
              sprintf(filenam,"/tmp/output_f0_%d.txt",num_file_);
	      fo=fopen(filenam,"w");
	      for (count=0;count<KEEP_SIZE;count++) 
                  fprintf(fo,"%ld %f %f\n",ch,bufout0_sta[count].real(),bufout0_sta[count].imag());
	      for (count=0;count<KEEP_SIZE;count++) 
                  fprintf(fo,"%ld %f %f\n",ch,bufout0_sto[count].real(),bufout0_sto[count].imag());
	      fclose(fo);
              sprintf(filenam,"/tmp/output_favg%ld_%d.txt",ch,num_file_);
	      fo=fopen(filenam,"w");
	      for (count=0;count<KEEP_SIZE*2;count++) 
                  fprintf(fo,"%ld %f %f\n",ch,spoofing_average[count].real(),spoofing_average[count].imag());
	      fclose(fo);
              num_file_++;
              printf("\n***********\n");
	      fo=fopen("/tmp/output_t1.txt","a");
	      for (count=0;count<CHUNK_SIZE;count++) fprintf(fo,"%ld %f %f\n",ch,in[count].real(),in[count].imag());
	      fclose(fo);
	      fo=fopen("/tmp/output_t2.txt","a");
	      for (count=0;count<CHUNK_SIZE;count++) fprintf(fo,"%ld %f %f\n",ch,carre[count].real(),carre[count].imag());
	      fclose(fo);
	     }
*/
        }
    if (avg_index_==Navg)    // restart averaging
      {
// https://www.libvolk.org/doxygen/volk_32fc_index_max_16u.html
// Finds and returns the index which contains the maximum magnitude for complex points in the given vector
       volk_32fc_index_max_16u(&maxpos, &spoofing_average_mul[KEEP_SIZE-10],20);  // max value where there should be no satellite
       maxvallim=10.*norm(spoofing_average_mul[KEEP_SIZE-10+maxpos]);   // abs(max)^2  TODO adjust x10 ?
       count=0;
       meanarg=0.;
       meanabs=0.;
#ifdef avgcpl
       meandiv={0.,0.};
#endif
       do
         {volk_32fc_index_max_16u(&maxpos, spoofing_average_mul, KEEP_SIZE*2); 
//          printf("%hd:\tangle=%.3f\t-\t",maxpos,arg(spoofing_average_mul[maxpos])); // atan2(spoofing_average_mul[maxpos].imag(),spoofing_average_mul[maxpos].real()));
//          printf("div: mag=%.3f angle=%.3f\n", abs(spoofing_average_div[maxpos]),arg(spoofing_average_div[maxpos]));
          maxval=norm(spoofing_average_mul[maxpos]);
          spoofing_average_mul[maxpos]=gr_complex{0.,0.}; //  on met ce bin a 0 et on itere sur ses voisins
          if (maxpos>0)
             spoofing_average_mul[maxpos-1]=gr_complex{0.,0.}; //  on met ce bin a 0 et on itere sur ses voisins
          if (maxpos>1)
             spoofing_average_mul[maxpos-2]=gr_complex{0.,0.}; //  on met ce bin a 0 et on itere sur ses voisins
          if (maxpos<2*KEEP_SIZE-1)
             spoofing_average_mul[maxpos+1]=gr_complex{0.,0.}; //  on met ce bin a 0 et on itere sur ses voisins
          if (maxpos<2*KEEP_SIZE-2)
             spoofing_average_mul[maxpos+2]=gr_complex{0.,0.}; //  on met ce bin a 0 et on itere sur ses voisins
#ifdef avgcpl
          meandiv+=spoofing_average_div[maxpos];
          stddiv[count]=spoofing_average_div[maxpos];
#else
          meanarg+=arg(spoofing_average_div[maxpos]);
          meanabs+=abs(spoofing_average_div[maxpos]);
          stdarg[count]=arg(spoofing_average_div[maxpos]);
          stdabs[count]=abs(spoofing_average_div[maxpos]);
#endif
          count++;
         }
       while ((maxval>maxvallim)&&(count<MAXSAT));
       if (count>0)
          {meanabs/=(float)(count);
           meanarg/=(float)(count);
#ifdef avgcpl
           meandiv/=(float)(count);
#endif
           stdargres_=0.;
#ifdef avgcpl
           for (i=0;i<count;i++) {stdargres_+=arg(stddiv[i]-meandiv)*arg(stddiv[i]-meandiv);}
#else
           for (i=0;i<count;i++) {stdargres_+=(stdarg[i]-meanarg)*(stdarg[i]-meanarg);}
#endif
           stdargres_/=(float)(count);
           //weight={-sqrtf(sqrtf(meanmod))*cosf(meanarg/2),-sqrtf(sqrtf(meanmod))*sinf(meanarg/2)}; // norm=|.|^2 and .=(A1/A2)^2  WHY +pi ?
           if (stdargres_<=STD_THRESHOLD) //  spoofing
             {weightabs=0.;
              weightarg=0.;
              c=0;
              cnt=0;
              printf("selected: ");
              do
                {if (abs(stdarg[c]-meanarg)<0.1) {printf("%d ",c);weightarg+=stdarg[c];weightabs+=stdabs[c];cnt++;} // only keep spoofing SV
                 c++;
                }
              while ((c<count)&&(c<MAXKEEP));
              weightarg/=(float)c;
              weightabs/=(float)c;
              printf("\n");
#ifdef avgcpl
              weight_={sqrtf(abs(meandiv))*cosf(arg(meandiv)/2),sqrtf(abs(meandiv))*sinf(arg(meandiv)/2)}; // .=(A1/A2)^2  WHY +pi ?
#else
              weight_={sqrtf(weightabs)*cosf(weightarg/2),sqrtf(weightabs)*sinf(weightarg/2)}; // .=(A1/A2)^2  WHY +pi ?
#endif
             }
           printf("%d:\tmeanarg=%.4f\tmeanabs=%.3f\tstdargres_=%.5f\tweightabs=%.2f,weightarg=%.2f",count,(meanarg),(meanabs),stdargres_,weightabs,weightarg);
           if (stdargres_>STD_THRESHOLD) printf("\n"); else printf(" /!\\\n");
          }
       memset(spoofing_average_mul,0,sizeof(gr_complex)*KEEP_SIZE*2); // start + stop
       memset(spoofing_average_div,0,sizeof(gr_complex)*KEEP_SIZE*2); // start + stop
       avg_index_=0;
      }
    // memcpy(output_items[0], input_items[0], noutput_items * input_signature()->sizeof_stream_item(ch));
    // memcpy(output_items[0], carre, noutput_items * input_signature()->sizeof_stream_item(ch));
    if (stdargres_>STD_THRESHOLD) // no spoofing
       {memcpy(output_items[0], input_items[0], noutput_items * input_signature()->sizeof_stream_item(ch));
       }
    else
      {volk_32fc_s32fc_multiply_32fc(carre,(const gr_complex*)input_items[0],weight_,CHUNK_SIZE); // -alpha*ant0
       volk_32fc_x2_add_32fc((gr_complex*)output_items[0], (const gr_complex*)input_items[1], carre, CHUNK_SIZE);
      }
//    delete plan;
    volk_free(carre);
    return noutput_items;
}
