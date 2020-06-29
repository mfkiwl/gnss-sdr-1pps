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

#include "jmf_spoofing.h"
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
#define STD_THRESHOLD 0.05 // rad
#define MAXSAT  20  // too many satellites will start detecting genuine constellation ?
#define MAXKEEP 7   // too many satellites will start detecting genuine constellation ?

#pragma message("JMF compile")

#define Navg 1  // FFT averages
#define moycpl  // average complex (if active) or average mag/phase (if inactive)

Gnss_Spoofing_Protect::Gnss_Spoofing_Protect(size_t sizeof_stream_item,
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue) : gr::sync_block("jmf_spoofing",
                               gr::io_signature::make(1, 20, sizeof_stream_item),
                               gr::io_signature::make(1, 1, sizeof_stream_item)),
                           d_ncopied_items(0),
                           d_queue(std::move(queue))
{
    printf("JMF Gnss_Spoofing_Protect\n");
/*
Ron Economos (April 5, 2020 10:58 AM)
To: discuss-gnuradio@gnu.org
I would use set_output_multiple() instead. See my previous e-mail for an 
example.
https://lists.gnu.org/archive/html/discuss-gnuradio/2019-08/msg00188.html
*/
    set_output_multiple(CHUNK_SIZE); // only trigger processing if that amount of samples was accumulated
    first_time_=0;
    memset(jmf_average_mul,0,sizeof(gr_complex)*KEEP_SIZE*2);
    memset(jmf_average_div,0,sizeof(gr_complex)*KEEP_SIZE*2);
    avg_index_=0;
    num_file_=0;
    weight_=1.;
    stdargres_=10.;
}


boost::shared_ptr<Gnss_Spoofing_Protect> gnss_sdr_make_spoof(size_t sizeof_stream_item, std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue)
{
//    unsigned int alignment = volk_get_alignment();
    boost::shared_ptr<Gnss_Spoofing_Protect> jmf_spoofing_(new Gnss_Spoofing_Protect(sizeof_stream_item, std::move(queue)));
    printf("JMF variable created\n");
/*
    jmf_average=(gr_complex*)volk_malloc(sizeof(gr_complex)*KEEP_SIZE*2,alignement); // 2* since sta followed by sto -> must fftshift to put 0 at center
    bufout0_sta=(gr_complex*)volk_malloc(sizeof(gr_complex)*KEEP_SIZE,alignement); 
    bufout0_sto=(gr_complex*)volk_malloc(sizeof(gr_complex)*KEEP_SIZE,alignement); 
    bufout_sta=(gr_complex*)volk_malloc(sizeof(gr_complex)*KEEP_SIZE,alignement); 
    bufout_sto=(gr_complex*)volk_malloc(sizeof(gr_complex)*KEEP_SIZE,alignement); 
*/
    return jmf_spoofing_;
}

int Gnss_Spoofing_Protect::work(int noutput_items,
    gr_vector_const_void_star &input_items,
    gr_vector_void_star &output_items)
{   long unsigned int ch = 0;
    int i,c,cnt;
    uint16_t maxpos;
    unsigned int alignment = volk_get_alignment();
    float maxval,maxvallim;
    gr_complex weightcpl={0.,0.};
    gr_complex* bufin;
#ifdef moycpl
    gr_complex stddiv[MAXSAT];
#else
     float meanarg,meanabs,stdarg[MAXSAT],stdabs[MAXSAT],weightabs=0.,weightarg=0.;
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
//    printf("JMF block: %d items, %ld out, %ld in\n",noutput_items,output_items.size(),input_items.size());
//           JMF block: 32768 items, 1 out, 2 in
    if (input_items.size()!=2) first_time_=1; // don't save if other than 2 input channels
/////////// CHECK FILES RUN AT THE SAME RATE
// select same file with different filenames
/*
    volk_32fc_s32fc_multiply_32fc(carre, (const gr_complex*)input_items[0],-1,   CHUNK_SIZE); 
    volk_32fc_x2_add_32fc        (carre, (const gr_complex*)input_items[1],carre,CHUNK_SIZE); 
    maxpos=0;
    for (i=0;i<CHUNK_SIZE;i++)
       {if (carre[i].real()!=0) {printf("%d: real !=0\n",i);maxpos=i;}
        if (carre[i].imag()!=0) {printf("%d: imag !=0\n",i);maxpos=i;}
       }   
if (maxpos!=0) {printf("JMF: sync error\n");fflush(stdout);}
*/
/////////// END CHECK FILES RUN AT THE SAME RATE
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
              volk_32fc_x2_add_32fc( jmf_average_mul,            jmf_average_mul,           bufout_sta,KEEP_SIZE); 
              volk_32fc_x2_add_32fc(&jmf_average_mul[KEEP_SIZE],&jmf_average_mul[KEEP_SIZE],bufout_sto,KEEP_SIZE); 

              volk_32fc_x2_divide_32fc(bufout_sta,bufout0_sta,plan->get_outbuf(),KEEP_SIZE);  // CH1/CH0
              volk_32fc_x2_divide_32fc(bufout_sto,bufout0_sto,&plan->get_outbuf()[CHUNK_SIZE-KEEP_SIZE-1],KEEP_SIZE); 
// for (i=0;i<KEEP_SIZE;i++) {bufout_sta[i]=sqrt(bufout_sta[i]); bufout_sto[i]=sqrt(bufout_sto[i]);}
              volk_32fc_x2_add_32fc( jmf_average_div,            jmf_average_div,           bufout_sta,KEEP_SIZE); 
              volk_32fc_x2_add_32fc(&jmf_average_div[KEEP_SIZE],&jmf_average_div[KEEP_SIZE],bufout_sto,KEEP_SIZE); 
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
                  fprintf(fo,"%ld %f %f\n",ch,jmf_average[count].real(),jmf_average[count].imag());
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
      {//volk_32fc_magnitude_squared_32f(jmf_mag,jmf_average,KEEP_SIZE*2);
// https://www.libvolk.org/doxygen/volk_32fc_index_max_16u.html
// Finds and returns the index which contains the maximum magnitude for complex points in the given vector
       volk_32fc_index_max_16u(&maxpos, &jmf_average_mul[KEEP_SIZE-10],20);  // max value where there should be no satellite
       maxvallim=10.*norm(jmf_average_mul[KEEP_SIZE-10+maxpos]);   // abs(max)^2  TODO adjust x10 ?
       count=0;
#ifdef moycpl
       meandiv={0.,0.};
       weightcpl={0.,0.};
#else
       meanarg=0.;
       meanabs=0.;
#endif
       do
         {volk_32fc_index_max_16u(&maxpos, jmf_average_mul, KEEP_SIZE*2); 
//          printf("%hd:\tangle=%.3f\t-\t",maxpos,arg(jmf_average_mul[maxpos])); // atan2(jmf_average_mul[maxpos].imag(),jmf_average_mul[maxpos].real()));
//          printf("div: mag=%.3f angle=%.3f\n", abs(jmf_average_div[maxpos]),arg(jmf_average_div[maxpos]));
          maxval=norm(jmf_average_mul[maxpos]);
          jmf_average_mul[maxpos]=gr_complex{0.,0.}; //  on met ce bin a 0 et on itere sur ses voisins
          if (maxpos>0)
             jmf_average_mul[maxpos-1]=gr_complex{0.,0.}; //  on met ce bin a 0 et on itere sur ses voisins
          if (maxpos>1)
             jmf_average_mul[maxpos-2]=gr_complex{0.,0.}; //  on met ce bin a 0 et on itere sur ses voisins
          if (maxpos<2*KEEP_SIZE-1)
             jmf_average_mul[maxpos+1]=gr_complex{0.,0.}; //  on met ce bin a 0 et on itere sur ses voisins
          if (maxpos<2*KEEP_SIZE-2)
             jmf_average_mul[maxpos+2]=gr_complex{0.,0.}; //  on met ce bin a 0 et on itere sur ses voisins
#ifdef moycpl
          meandiv+=sqrt(jmf_average_div[maxpos]); // racine des coefficient mis au carr'e
          stddiv[count]=sqrt(jmf_average_div[maxpos]);
#else
          meanarg+=(arg(jmf_average_div[maxpos])/2);
          meanabs+=sqrt(abs(jmf_average_div[maxpos]));
          stdarg[count]=(arg(jmf_average_div[maxpos])/2);
          stdabs[count]=sqrt(abs(jmf_average_div[maxpos]));
#endif
          count++;
         }
       while ((maxval>maxvallim)&&(count<MAXSAT));
       if (count>0)
          {
#ifdef moycpl
           meandiv/=(float)(count);
#else
           meanabs/=(float)(count);
           meanarg/=(float)(count);
#endif
           stdargres_=0.;
#ifdef moycpl
           for (i=0;i<count;i++) {stdargres_+=(arg(stddiv[i])-arg(meandiv))*(arg(stddiv[i])-arg(meandiv));} // ecart type sur les phases
//           for (i=0;i<count;i++) {printf("%d: %f\n",i,abs(stddiv[i]));}  // affiche le module de tous les poids
#else
           for (i=0;i<count;i++) {stdargres_+=(stdarg[i]-meanarg)*(stdarg[i]-meanarg);}
#endif
           stdargres_/=(float)(count);
           printf("%d:\t%.5f\t",count,stdargres_);fflush(stdout);
           if (stdargres_<=STD_THRESHOLD) //  spoofing
             {
#ifndef moycpl
              weightabs=0.;
              weightarg=0.;
#endif
              c=0;
              cnt=0;
              do {
#ifdef moycpl
                if ((arg(stddiv[c])-arg(meandiv))<0.1) 
                    {weightcpl+=stddiv[c];
                     cnt++;
                    }
#else
                if (abs(stdarg[c]-meanarg)<0.1) 
                    {printf("%d ",c); 
                     weightarg+=stdarg[c];
                     weightabs+=stdabs[c];
                     cnt++;
                    } // only keep spoofing SV
#endif
                c++;
              } while ((c<count)&&(cnt<MAXKEEP));
#ifdef moycpl
              weightcpl/=(float)cnt;
#else
              weightarg/=(float)cnt;
              weightabs/=(float)cnt;
#endif
#ifdef moycpl  // should be -sqrt, -sqrt to subtract, but with \pi phase rotation remove the '-'
              weight_=weightcpl; // .=(A1/A2)^2  WHY +pi ?
#else
              weight_={(weightabs)*cosf(weightarg),(weightabs)*sinf(weightarg)}; // .=(A1/A2)^2  WHY +pi ?
#endif
             }
           printf("weightabs=%.2f,weightarg=%.2f",abs(weight_),arg(weight_));
           if (stdargres_>STD_THRESHOLD) printf("\n"); else printf(" /!\\\n");
          }
       memset(jmf_average_mul,0,sizeof(gr_complex)*KEEP_SIZE*2); // start + stop
       memset(jmf_average_div,0,sizeof(gr_complex)*KEEP_SIZE*2); // start + stop
       avg_index_=0;
      }
    // memcpy(output_items[0], input_items[0], noutput_items * input_signature()->sizeof_stream_item(ch));
    // memcpy(output_items[0], carre, noutput_items * input_signature()->sizeof_stream_item(ch));
    if (stdargres_>STD_THRESHOLD) // no spoofing
       {memcpy(output_items[0], input_items[0], noutput_items * input_signature()->sizeof_stream_item(ch));
       }
    else
      {volk_32fc_s32fc_multiply_32fc(carre,(const gr_complex*)input_items[1],weight_,CHUNK_SIZE); // -alpha*ant0
       volk_32fc_x2_add_32fc((gr_complex*)output_items[0], (const gr_complex*)input_items[0], carre, CHUNK_SIZE);
      }
//    delete plan;
    volk_free(carre);
    return noutput_items;
}
