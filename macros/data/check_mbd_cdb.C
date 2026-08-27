
#ifndef _CHECK_MBD_CDB_C_
#define _CHECK_MBD_CDB_C_

#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libffamodules.so )

// Diagnostic for the "MinBiasCut rejects every event" runs: queries which
// MBD calibration CDB payloads resolve to a URL vs. come back empty
// ("MISSING") for a given run's timestamp -- no event processing needed,
// CDBInterface::getUrl() is a pure IOV lookup keyed on recoConsts
// TIMESTAMP, same mechanism Fun4All_Dijets_AuAu.C uses at line ~98.
// Pass neighbor_run (a run known to work, e.g. the nearest good run in
// the GRL) to print its resolution side by side for comparison.
void check_mbd_cdb(
    const int run,
    const int neighbor_run = -1,
    const std::string & cdbtag = "newcdbtag"
)
{
    static const std::vector<std::string> domains = {
        "MBD_QFIT", "MBD_SAMPMAX", "MBD_SLEWCORR", "MBD_T0CORR",
        "MBD_TIMECORR", "MBD_TQ_T0", "MBD_TT_T0"
    };

    auto * rc = recoConsts::instance();
    rc -> set_StringFlag( "CDB_GLOBALTAG", cdbtag );

    auto check = [&]( const int ts )
    {
        rc -> set_uint64Flag( "TIMESTAMP", ts );
        for ( const auto & domain : domains )
        {
            const std::string url = CDBInterface::instance() -> getUrl( domain );
            std::cout << "  run " << ts << "  " << domain << ": "
                       << ( url.empty() ? "MISSING" : url ) << std::endl;
        }
    };

    std::cout << "=== run " << run << " ===" << std::endl;
    check( run );

    if ( neighbor_run > 0 )
    {
        std::cout << "=== neighbor run " << neighbor_run << " (for comparison) ===" << std::endl;
        check( neighbor_run );
    }
}

#endif
