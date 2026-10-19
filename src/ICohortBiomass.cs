using Landis.Core;
using Landis.Library.UniversalCohorts;
using Landis.SpatialModeling;
using Landis.Utilities;

namespace Landis.Extension.Succession.BiomassPnET
{
    public interface ICohortBiomass : ICalculator
    {
        /// <summary>
        /// Computes the change in an individual cohort's biomass due to annual
        /// growth and mortality.
        /// </summary>
        int ComputeChange(ICohort cohort, ActiveSite site);

        /// <summary>
        /// Computes the percentage of a cohort's biomass that is non-woody.
        /// </summary>
        Percentage ComputeNonWoodyPercentage(ICohort cohort, ActiveSite site);
    }
}
