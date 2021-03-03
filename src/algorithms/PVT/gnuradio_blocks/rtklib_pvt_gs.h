/*!
 * \file rtklib_pvt_gs.h
 * \brief Interface of a Position Velocity and Time computation block
 * \author Javier Arribas, 2017. jarribas(at)cttc.es
 *
 * -----------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2020  (see AUTHORS file for a list of contributors)
 *
 * GNSS-SDR is a software defined Global Navigation
 *          Satellite Systems receiver
 *
 * This file is part of GNSS-SDR.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * -----------------------------------------------------------------------------
 */

#ifndef GNSS_SDR_RTKLIB_PVT_GS_H
#define GNSS_SDR_RTKLIB_PVT_GS_H

#include "gnss_synchro.h"
#include "rtklib.h"
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <gnuradio/sync_block.h>  // for sync_block
#include <gnuradio/types.h>       // for gr_vector_const_void_star
#include <pmt/pmt.h>              // for pmt_t
#include <chrono>                 // for system_clock
#include <cstddef>                // for size_t
#include <cstdint>                // for int32_t
#include <ctime>                  // for time_t
#include <map>                    // for map
#include <memory>                 // for shared_ptr, unique_ptr
#include <string>                 // for string
#include <sys/types.h>            // for key_t
#include <vector>                 // for vector
#if GNURADIO_USES_STD_POINTERS
#else
#include <boost/shared_ptr.hpp>
#endif

#include "vxi11_user.h"

class Beidou_Dnav_Almanac;
class Beidou_Dnav_Ephemeris;
class Galileo_Almanac;
class Galileo_Ephemeris;
class GeoJSON_Printer;
class Gps_Almanac;
class Gps_Ephemeris;
class Gpx_Printer;
class Kml_Printer;
class Monitor_Pvt_Udp_Sink;
class Nmea_Printer;
class Pvt_Conf;
class Rinex_Printer;
class Rtcm_Printer;
class Rtklib_Solver;
class rtklib_pvt_gs;

#if GNURADIO_USES_STD_POINTERS
using rtklib_pvt_gs_sptr = std::shared_ptr<rtklib_pvt_gs>;
#else
using rtklib_pvt_gs_sptr = boost::shared_ptr<rtklib_pvt_gs>;
#endif

rtklib_pvt_gs_sptr rtklib_make_pvt_gs(uint32_t nchannels,
    const Pvt_Conf& conf_,
    const rtk_t& rtk,
    const double PPS_Kp,
    const double PPS_Ki, 
    const bool SMA_internal_source_clock, 
    double LO_external_frequ,
    const bool PPS_correction,
    const bool PPS_estimator_selected,
    const std::string SMA_IP_address);

/*!
 * \brief This class implements a block that computes the PVT solution using the RTKLIB integrated library
 */
class rtklib_pvt_gs : public gr::sync_block
{
public:
    ~rtklib_pvt_gs();  //!< Default destructor

    /*!
     * \brief Get latest set of GPS ephemeris from PVT block
     */
    std::map<int, Gps_Ephemeris> get_gps_ephemeris_map() const;

    /*!
     * \brief Get latest set of GPS almanac from PVT block
     */
    std::map<int, Gps_Almanac> get_gps_almanac_map() const;

    /*!
     * \brief Get latest set of Galileo ephemeris from PVT block
     */
    std::map<int, Galileo_Ephemeris> get_galileo_ephemeris_map() const;

    /*!
     * \brief Get latest set of Galileo almanac from PVT block
     */
    std::map<int, Galileo_Almanac> get_galileo_almanac_map() const;

    /*!
     * \brief Get latest set of BeiDou DNAV ephemeris from PVT block
     */
    std::map<int, Beidou_Dnav_Ephemeris> get_beidou_dnav_ephemeris_map() const;

    /*!
     * \brief Get latest set of BeiDou DNAV almanac from PVT block
     */
    std::map<int, Beidou_Dnav_Almanac> get_beidou_dnav_almanac_map() const;

    /*!
     * \brief Clear all ephemeris information and the almanacs for GPS and Galileo
     */
    void clear_ephemeris();

    /*!
     * \brief Get the latest Position WGS84 [deg], Ground Velocity, Course over Ground, and UTC Time, if available
     */
    bool get_latest_PVT(double* longitude_deg,
        double* latitude_deg,
        double* height_m,
        double* ground_speed_kmh,
        double* course_over_ground_deg,
        time_t* UTC_time) const;

    int work(int noutput_items, gr_vector_const_void_star& input_items,
        gr_vector_void_star& output_items);  //!< PVT Signal Processing

private:
    friend rtklib_pvt_gs_sptr rtklib_make_pvt_gs(uint32_t nchannels,
        const Pvt_Conf& conf_,
        const rtk_t& rtk,
	const double PPS_Kp,
	const double PPS_Ki,
	const bool SMA_internal_source_clock,
	double LO_external_frequ,
	const bool PPS_correction,
	const bool PPS_estimator_selected,
	const std::string SMA_IP_address);

    rtklib_pvt_gs(uint32_t nchannels,
        const Pvt_Conf& conf_,
        const rtk_t& rtk,
	const double PPS_Kp,
	const double PPS_Ki,
	const bool SMA_internal_source_clock,
	double LO_external_frequ,
	const bool PPS_correction,
	const bool PPS_estimator_selected,
	const std::string SMA_IP_address);


    void msg_handler_telemetry(const pmt::pmt_t& msg);

    void initialize_and_apply_carrier_phase_offset();

    void apply_rx_clock_offset(std::map<int, Gnss_Synchro>& observables_map,
        double rx_clock_offset_s);

    std::map<int, Gnss_Synchro> interpolate_observables(const std::map<int, Gnss_Synchro>& observables_map_t0,
        const std::map<int, Gnss_Synchro>& observables_map_t1,
        double rx_time_s);

    inline std::time_t convert_to_time_t(const boost::posix_time::ptime pt) const
    {
        return (pt - boost::posix_time::ptime(boost::gregorian::date(1970, 1, 1))).total_seconds();
    }

    std::vector<std::string> split_string(const std::string& s, char delim) const;

    typedef struct
    {
        long mtype;  // NOLINT(google-runtime-int) required by SysV queue messaging
        double ttff;
    } d_ttff_msgbuf;
    bool send_sys_v_ttff_msg(d_ttff_msgbuf ttff);

    bool save_gnss_synchro_map_xml(const std::string& file_name);  // debug helper function
    bool load_gnss_synchro_map_xml(const std::string& file_name);  // debug helper function

    std::shared_ptr<Rtklib_Solver> d_internal_pvt_solver;
    std::shared_ptr<Rtklib_Solver> d_user_pvt_solver;

    std::unique_ptr<Rinex_Printer> d_rp;
    std::unique_ptr<Kml_Printer> d_kml_dump;
    std::unique_ptr<Gpx_Printer> d_gpx_dump;
    std::unique_ptr<Nmea_Printer> d_nmea_printer;
    std::unique_ptr<GeoJSON_Printer> d_geojson_printer;
    std::unique_ptr<Rtcm_Printer> d_rtcm_printer;
    std::unique_ptr<Monitor_Pvt_Udp_Sink> d_udp_sink_ptr;

    std::chrono::time_point<std::chrono::system_clock> d_start;
    std::chrono::time_point<std::chrono::system_clock> d_end;

    std::string d_dump_filename;
    std::string d_xml_base_path;
    std::string d_local_time_str;

    std::vector<bool> d_channel_initialized;
    std::vector<double> d_initial_carrier_phase_offset_estimation_rads;

    enum StringValue_
    {
        evGPS_1C,
        evGPS_2S,
        evGPS_L5,
        evSBAS_1C,
        evGAL_1B,
        evGAL_5X,
        evGLO_1G,
        evGLO_2G,
        evBDS_B1,
        evBDS_B2,
        evBDS_B3
    };
    std::map<std::string, StringValue_> d_mapStringValues;
    std::map<int, Gnss_Synchro> d_gnss_observables_map;
    std::map<int, Gnss_Synchro> d_gnss_observables_map_t0;
    std::map<int, Gnss_Synchro> d_gnss_observables_map_t1;

    boost::posix_time::time_duration d_utc_diff_time;

    size_t d_gps_ephemeris_sptr_type_hash_code;
    size_t d_gps_iono_sptr_type_hash_code;
    size_t d_gps_utc_model_sptr_type_hash_code;
    size_t d_gps_cnav_ephemeris_sptr_type_hash_code;
    size_t d_gps_cnav_iono_sptr_type_hash_code;
    size_t d_gps_cnav_utc_model_sptr_type_hash_code;
    size_t d_gps_almanac_sptr_type_hash_code;
    size_t d_galileo_ephemeris_sptr_type_hash_code;
    size_t d_galileo_iono_sptr_type_hash_code;
    size_t d_galileo_utc_model_sptr_type_hash_code;
    size_t d_galileo_almanac_helper_sptr_type_hash_code;
    size_t d_galileo_almanac_sptr_type_hash_code;
    size_t d_glonass_gnav_ephemeris_sptr_type_hash_code;
    size_t d_glonass_gnav_utc_model_sptr_type_hash_code;
    size_t d_glonass_gnav_almanac_sptr_type_hash_code;
    size_t d_beidou_dnav_ephemeris_sptr_type_hash_code;
    size_t d_beidou_dnav_iono_sptr_type_hash_code;
    size_t d_beidou_dnav_utc_model_sptr_type_hash_code;
    size_t d_beidou_dnav_almanac_sptr_type_hash_code;

    double d_rinex_version;
    double d_rx_time;

    key_t d_sysv_msg_key;
    int d_sysv_msqid;

    int32_t d_rinexobs_rate_ms;
    int32_t d_rtcm_MT1045_rate_ms;  // Galileo Broadcast Ephemeris
    int32_t d_rtcm_MT1019_rate_ms;  // GPS Broadcast Ephemeris (orbits)
    int32_t d_rtcm_MT1020_rate_ms;  // GLONASS Broadcast Ephemeris (orbits)
    int32_t d_rtcm_MT1077_rate_ms;  // The type 7 Multiple Signal Message format for the USA’s GPS system, popular
    int32_t d_rtcm_MT1087_rate_ms;  // GLONASS MSM7. The type 7 Multiple Signal Message format for the Russian GLONASS system
    int32_t d_rtcm_MT1097_rate_ms;  // Galileo MSM7. The type 7 Multiple Signal Message format for Europe’s Galileo system
    int32_t d_rtcm_MSM_rate_ms;
    int32_t d_kml_rate_ms;
    int32_t d_gpx_rate_ms;
    int32_t d_geojson_rate_ms;
    int32_t d_nmea_rate_ms;
    int32_t d_last_status_print_seg;  // for status printer
    int32_t d_output_rate_ms;
    int32_t d_display_rate_ms;
    int32_t d_report_rate_ms;
    int32_t d_max_obs_block_rx_clock_offset_ms;

    uint32_t d_nchannels;
    uint32_t d_type_of_rx;

    bool d_dump;
    bool d_dump_mat;
    bool d_rinex_output_enabled;
    bool d_rinex_header_written;
    bool d_rinex_header_updated;
    bool d_geojson_output_enabled;
    bool d_gpx_output_enabled;
    bool d_kml_output_enabled;
    bool d_nmea_output_file_enabled;
    bool d_first_fix;
    bool d_xml_storage;
    bool d_flag_monitor_pvt_enabled;
    bool d_show_local_time_zone;
    bool d_waiting_obs_block_rx_clock_offset_correction_msg;
    bool d_enable_rx_clock_correction;
    bool d_rtcm_writing_started;
    bool d_rtcm_enabled;
  
    bool d_PPS_correction; 
    bool d_PPS_estimator_selected;
    double d_PPS_Kp;
    double d_PPS_Ki;
    double d_LO_external_frequ;
    double d_LO_external_frequ_init;

    std::string d_SMA_IP_address;
    VXI11_CLINK *d_clink;


//-------------1PPS variables 

    double d_FREQU_correction;
    double d_pps_prev_error;
    double d_estimator;
    double d_estimator_prev;
    double d_pps_init_offset;
    double d_pps_offset;
};

#endif  // GNSS_SDR_RTKLIB_PVT_GS_H
