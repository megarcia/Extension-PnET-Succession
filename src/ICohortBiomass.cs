using Landis.Core;
using Landis.Library.UniversalCohorts;
using Landis.SpatialModeling;
using Landis.Utilities;

namespace Landis.Extension.Succession.BiomassPnET
{
    public interface ICohortBiomass : Library.UniversalCohorts.ICalculator
    {
        /// <summary>
        /// Computes the change in an individual cohort's biomass due to annual
        /// growth and mortality.
        /// </summary>
        int ComputeChange(Library.UniversalCohorts.ICohort cohort, ActiveSite site);

        /// <summary>
        /// Computes the non-woody percentage of a cohort's biomass.
        /// </summary>
        Percentage ComputeNonWoodyPercentage(Library.UniversalCohorts.ICohort cohort, ActiveSite site);
    }
}
