﻿/// the only class and method in this file are unreferenced -- 
/// will be deprecated and removed.

using Landis.Core;
using Landis.Library.PnETCohorts;
using Landis.SpatialModeling;

namespace Landis.Extension.Succession.BiomassPnET.DisturbanceReductions
{
    class Events
    {
        //---------------------------------------------------------------------
        public static void CohortDied(object sender, Landis.Library.UniversalCohorts.DeathEventArgs eventArgs)
        {
            ExtensionType disturbanceType = eventArgs.DisturbanceType;
        }
        //---------------------------------------------------------------------

    }
}
