// Matthew Garcia: (10.09.2025)
//
// There's no reason for this method to have its own namespace, 
// complicating the PlugIn code. I have renamed this class and file
// "Mortality," since I wanted to establish such a class anyway 
// (especially for use with related methods in the PnET-Cohort 
// Library), and reverted it to the Succession.BiomassPnET 
// namespace.

// NOTE: DeathEventArgs --> Library.UniversalCohorts 
// NOTE: ExtensionType --> Landis.Core

using Landis.Core;
using Landis.Library.PnETCohorts;
using Landis.Library.UniversalCohorts;
// using Landis.SpatialModeling;

namespace Landis.Extension.Succession.BiomassPnET
{
    class Mortality
    {
        public static void CohortDied(object sender, DeathEventArgs eventArgs)
        {
            ExtensionType disturbanceType = eventArgs.DisturbanceType;
        }
    }
}
