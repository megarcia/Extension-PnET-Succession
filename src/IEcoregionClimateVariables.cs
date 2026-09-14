//MG20260914 This interface was part of the original way to feed climate
//           information to PnET-Succession but is now redundant with the
//           new LANDIS-II Climate Library and causes several redundancies
//           in the PnET Cohort Library. This interface and its counterpart
//           functionality in the PnET Cohort Library will be deprecated. 

using Landis.Library.PnETCohorts;

namespace Landis.Extension.Succession.BiomassPnET
{
    public interface IEcoregionClimateVariables
    {

        float PAR0 { get; }
        float Prec { get; }
        float Tday { get; }
        float VPD { get; }
        float Year { get; }
        float DaySpan { get; }
        float Daylength { get; }
        byte Month { get; }
        float Tave { get; }
        float Tmin { get; }
        float Tmax { get; }
        float CO2 { get; }
        float O3 { get; }
        float SPEI { get; }
        
        SpeciesPnETVariables this[string species] { get; }
    }
}
