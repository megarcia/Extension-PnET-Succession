// Original Authors: Arjan de Bruijn, Brian R. Miranda

// John McNabb: (02.04.2019)
//
// Summary of changes to allow the climate library to be used with PnET-Succession:
//   (1) Added ClimateRegionData class based on that of NECN to hold the climate library data. This is Initialized by a call
//       to InitialClimateLibrary() in Plugin.Initialize().
//   (2) Modified EcoregionPnET to add GetClimateRegionData() which grabs climate data from ClimateRegionData.  This uses an intermediate
//       MonthlyClimateRecord instance which is similar to ObservedClimate.
//       MonthlyClimateRecord instance which is similar to ObservedClimate.
//   (3) Added ClimateRegionPnETVariables class which is a copy of the EcoregionPnETVariables class which uses MonthlyClimateRecord rather than
//       ObservedClimate. I had hoped to use the same class, but the definition of IObservedClimate prevents MonthlyClimateRecord from implementing it.
//       IMPORTANT NOTE: The climate library precipation is in cm/month, so that it is converted to mm/month in MonthlyClimateRecord.
//   (4) Modified Plugin.AgeCohorts() and SiteCohorts.SiteCohorts() to call either EcoregionPnET.GetClimateRegoinData() or EcoregionPnET.GetData()
//       depending on whether the climate library is enabled.
//
// Enabling the climate library with PnET:
//   (1) Indicate the climate library configuration file in the 'PnET-succession' configuration file using the 'ClimateConfigFile' parameter, e.g.
//        ClimateConfigFile	"./climate-generator-baseline.txt"
//
// NOTE: Use of the climate library is OPTIONAL.  If the 'ClimateConfigFile' parameter is missing (or commented-out) of the 'PnET-succession'
// configuration file, then PnET reverts to using climate data as specified by the 'ClimateFileName' column in the 'EcoregionParameters' file
// given in the 'PnET-succession' configuration file.
//
// NOTE: This uses a version (v4?) of the climate library that exposes AnnualClimate_Monthly.MonthlyOzone[] and .MonthlyCO2[].

////////

// Matthew Garcia: (10.09.2026)
//
// Many of John's comments above contain file and class names that I
// have modified in recent editing of the PnET-Cohort Library. I've
// also added numerous hints here on where to find many of the classes
// invoked in this PlugIn, since they're all over the place in Landis.
//
// The PnET-Cohorts Library is used only by the PnET-Succession 
// Extension and not by any other Extension in Landis. Therefore,
// I plan to update this PlugIn first so that it's consistent with
// the updated PnET-Cohorts Library (Extension v6.1), then actually
// merge the PnET-Cohorts Library with the PnET-Succession Extension 
// so that they're one package within Landis (Extension v6.2). There
// is a lot of duplication and repetition from PnET-Cohorts in the 
// Extension/PlugIn code that can be eliminated, and other calls and 
// procedures can be streamlined somewhat with this merger. It will 
// also reduce the compile procedure from two separate steps (first
// the PnET - Cohorts Library, then the actual Extension) to one 
// compiler step.
//
// NOTE: ActiveSite --> Landis.SpatialModeling
// NOTE: Climate --> Landis.Library.Climate
// NOTE: Cohort --> Library.PnETCohorts
// NOTE: Constants --> Library.PnETCohorts
// NOTE: DeathEventArgs --> Library.UniversalCohorts 
// NOTE: Directory --> Landis.Utilities
// NOTE: ExtensionType --> Landis.Core
// NOTE: Globals --> Library.PnETCohorts
// NOTE: Hydrology --> Library.PnETCohorts
// NOTE: Hydrology_SaxtonRawls --> Library.PnETCohorts
// NOTE: ICore --> Landis.Core
// NOTE: IDataset --> Landis.Library.InitialCommunities.Universal
// NOTE: IInputRaster --> Landis.SpatialModeling
// NOTE: IPnETEcoregionData --> Library.PnETCohorts
// NOTE: IHydrology --> Library.PnETCohorts
// NOTE: ISiteVar --> Landis.SpatialModeling
// NOTE: Names --> Library.PnETCohorts
// NOTE: ObservedClimate --> Library.PnETCohorts
// NOTE: PnETEcoregionData --> Library.PnETCohorts
// NOTE: PnETSpecies --> Library.PnETCohorts
// NOTE: Reproduction --> Library.Succession
// NOTE: SeedingAlgorithms --> Library.Succession
// NOTE: SiteCohorts --> Library.PnETCohorts
// NOTE: SiteVars --> Library.PnETCohorts
// NOTE: Soils --> Library.PnETCohorts
// NOTE: SpeciesParameters --> Library.PnETCohorts

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Diagnostics;
using Landis.Core;
using Landis.Library.Climate;
using Landis.Library.InitialCommunities.Universal;
using Landis.Library.PnETCohorts;
using Landis.Library.Succession;
using Landis.Library.Succession.DensitySeeding;
using Landis.SpatialModeling;
using Landis.Library.UniversalCohorts;
using Landis.Utilities;

namespace Landis.Extension.Succession.BiomassPnET
{
    public class PlugIn : Landis.Library.Succession.ExtensionBase 
    {
        public static PnETSpecies PnETSpecies;
        public static DateTime Date;
        public static ICore ModelCore;
        private static DateTime StartDate;
        private static Dictionary<ActiveSite, string> SiteOutputNames;
        public static float FTimeStep;
        public static bool UsingClimateLibrary;
        private Dictionary<ActiveSite, ICommunity> sitesAndCommunities;
        public static string InitialCommunitiesSpinup;
        public static int CohortBinSize;
        public static int ParallelThreads;
        private static readonly object threadLock = new object();
        private Dictionary<ActiveSite, uint> allKeys;
        public static float MinFolRatioFactor;
        MyClock m = null;

        public void DeathEvent(object sender, DeathEventArgs eventArgs)
        {
            ExtensionType disturbanceType = eventArgs.DisturbanceType;
            if (disturbanceType != null)
            {
                ActiveSite site = eventArgs.Site;
                if (disturbanceType.IsMemberOf("disturbance:fire"))
                    Reproduction.CheckForPostFireRegen(eventArgs.Cohort, site);
                else
                    Reproduction.CheckForResprouting(eventArgs.Cohort, site);
            }
        }

        string PnETDefaultsFolder
        {
            get
            {
                string defaultPath = Path.Combine(Path.GetDirectoryName(Reflection.Assembly.GetExecutingAssembly().Location), "Defaults");
                // If Linux, correct the path string
                if (Runtime.InteropServices.RuntimeInformation.IsOSPlatform(Runtime.InteropServices.OSPlatform.Linux))
                    defaultPath = defaultPath.Replace('\\', '/');
                return defaultPath;
            }
        }

        public PlugIn() : base(Names.ExtensionName)
        {
            LocalOutput.PnETOutputSites = Names.PnETOutputSites;
            // The number of thread workers to use in succession routines that have been optimized. Should
            // more or less match the number of cores in the computer thats running LANDIS-II's processor
            // this.ThreadCount = 3;
            // this.ThreadCount = 1;
            allKeys = new Dictionary<ActiveSite, uint>();
            sitesAndCommunities = new Dictionary<ActiveSite, ICommunity>();
        }

        public override void LoadParameters(string InputParameterFile, ICore mCore)
        {
            ModelCore = mCore;
            Names.parameters.Add(Names.ExtensionName, new Parameter<string>(Names.ExtensionName, InputParameterFile));
            //-------------PnET-Succession input files
            Dictionary<string, Parameter<string>> InputParameters = Names.LoadTable(Names.ExtensionName, Names.AllNames, null, true);
            InputParameters.ToList().ForEach(x => Names.parameters.Add(x.Key, x.Value));
            //-------------Read Species parameters input file
            List<string> SpeciesNames = ModelCore.Species.ToList().Select(x => x.Name).ToList();
            List<string> SpeciesParamNames = PnETSpecies.ParameterNames;
            SpeciesParamNames.Add(Names.PnETSpeciesParameters);
            Dictionary<string, Parameter<string>> speciesparameters = Names.LoadTable(Names.PnETSpeciesParameters, SpeciesNames, SpeciesParamNames);
            foreach (string key in speciesparameters.Keys)
                if (Names.parameters.ContainsKey(key))
                    throw new Exception("Parameter " + key + " was provided twice");
            speciesparameters.ToList().ForEach(x => Names.parameters.Add(x.Key, x.Value));
            //-------------Ecoregion parameters
            List<string> EcoregionNames = ModelCore.Ecoregions.ToList().Select(x => x.Name).ToList();
            List<string> EcoregionParameters = PnETEcoregionData.ParameterNames;
            Dictionary<string, Parameter<string>> ecoregionparameters = Names.LoadTable(Names.EcoregionParameters, EcoregionNames, EcoregionParameters);
            foreach (string key in ecoregionparameters.Keys)
                if (Names.parameters.ContainsKey(key))
                    throw new Exception("Parameter " + key + " was provided twice");
            ecoregionparameters.ToList().ForEach(x => Names.parameters.Add(x.Key, x.Value));
            //-------------DisturbanceReductionsParameterFile
            Parameter<string> DisturbanceReductionsParameterFile;
            if (Names.TryGetParameter(Names.DisturbanceReductions, out DisturbanceReductionsParameterFile))
            {
                Allocation.Initialize(DisturbanceReductionsParameterFile.Value, Names.parameters);
                Cohort.AgeOnlyDeathEvent += Mortality.CohortDied;
            }
            //---------------SaxtonAndRawlsParameterFile
            if (Names.parameters.ContainsKey(PressureHeadSaxton_Rawls.SaxtonAndRawlsParameters) == false)
            {
                Parameter<string> SaxtonAndRawlsParameterFile = new Parameter<string>(PressureHeadSaxton_Rawls.SaxtonAndRawlsParameters, (string)PnETDefaultsFolder + System.IO.Path.DirectorySeparatorChar + "SaxtonAndRawlsParameters.txt");
                Names.parameters.Add(PressureHeadSaxton_Rawls.SaxtonAndRawlsParameters, SaxtonAndRawlsParameterFile);
            }
            Dictionary<string, Parameter<string>> SaxtonAndRawlsParameters = Names.LoadTable(PressureHeadSaxton_Rawls.SaxtonAndRawlsParameters, null, PressureHeadSaxton_Rawls.ParameterNames);
            foreach (string key in SaxtonAndRawlsParameters.Keys)
                if (Names.parameters.ContainsKey(key))
                    throw new Exception("Parameter " + key + " was provided twice");
            SaxtonAndRawlsParameters.ToList().ForEach(x => Names.parameters.Add(x.Key, x.Value));
            //--------------PnETGenericParameterFile
            //----------See if user supplied overwriting default parameters
            List<string> RowLabels = new List<string>(Names.AllNames);
            RowLabels.AddRange(PnETSpecies.ParameterNames);
            if (Names.parameters.ContainsKey(Names.PnETGenericParameters))
            {
                Dictionary<string, Parameter<string>> genericparameters = Names.LoadTable(Names.PnETGenericParameters, RowLabels, null, true);
                foreach (KeyValuePair<string, Parameter<string>> par in genericparameters)
                {
                    if (Names.parameters.ContainsKey(par.Key)) throw new System.Exception("Parameter " + par.Key + " was provided twice");
                    Names.parameters.Add(par.Key, par.Value);
                }
            }
            //----------Load in default parameters to fill the gaps
            Parameter<string> PnETGenericDefaultParameterFile = new Parameter<string>(Names.PnETGenericDefaultParameters, (string)PnETDefaultsFolder + System.IO.Path.DirectorySeparatorChar + "PnETGenericDefaultParameters.txt");
            Names.parameters.Add(Names.PnETGenericDefaultParameters, PnETGenericDefaultParameterFile);
            Dictionary<string, Parameter<string>> genericdefaultparameters = Names.LoadTable(Names.PnETGenericDefaultParameters, RowLabels, null, true);
            foreach (KeyValuePair<string, Parameter<string>> par in genericdefaultparameters)
            {
                if (Names.parameters.ContainsKey(par.Key) == false)
                    Names.parameters.Add(par.Key, par.Value);
            }
            SiteOutputNames = new Dictionary<ActiveSite, string>();
            Parameter<string> OutputSitesFile;
            if (Names.TryGetParameter(LocalOutput.PnETOutputSites, out OutputSitesFile))
            {
                Dictionary<string, Parameter<string>> outputfiles = Names.LoadTable(LocalOutput.PnETOutputSites, null, AssignOutputFiles.ParameterNames.AllNames, true);
                AssignOutputFiles.MapCells(outputfiles, ref SiteOutputNames);
            }
        }

        public override void Initialize()
        {
            ModelCore.UI.WriteLine("Initializing " + Names.ExtensionName + " version " + typeof(PlugIn).Assembly.GetName().Version);
            Cohort.DeathEvent += DeathEvent;
            StartDate = new DateTime(((Parameter<int>)Names.GetParameter(Names.StartYear)).Value, 1, 15);
            Globals.InitializeCore(ModelCore, ((Parameter<ushort>)Names.GetParameter(Names.IMAX)).Value, StartDate);
            PnETPnETEcoregionData.Initialize();
            SiteVars.Initialize();
            Directory.EnsureExists("output");
            Timestep = ((Parameter<int>)Names.GetParameter(Names.Timestep)).Value;
            Parameter<string> CohortBinSizeParm = null;
            if (Names.TryGetParameter(Names.CohortBinSize, out CohortBinSizeParm))
            {
                if (Int32.TryParse(CohortBinSizeParm.Value, out CohortBinSize))
                {
                    if (CohortBinSize < Timestep)
                        throw new System.Exception("CohortBinSize cannot be smaller than Timestep.");
                    else
                        ModelCore.UI.WriteLine("  Succession timestep = " + Timestep + "; CohortBinSize = " + CohortBinSize + ".");
                }
                else
                    throw new System.Exception("CohortBinSize is not an integer value.");
            }
            else
                CohortBinSize = Timestep;
            string Parallel = ((Parameter<string>)Names.GetParameter(Names.Parallel)).Value;
            if (Parallel == "false")
            {
                ParallelThreads = 1;
                ModelCore.UI.WriteLine("  MaxParallelThreads = " + ParallelThreads.ToString() + ".");
            }
            else if (Parallel == "true")
            {
                ParallelThreads = -1;
                ModelCore.UI.WriteLine("  MaxParallelThreads determined by system.");
            }
            else
            {
                if (Int32.TryParse(Parallel, out ParallelThreads))
                {
                    if (ParallelThreads < 1)
                        throw new System.Exception("Parallel cannot be < 1.");
                    else
                        ModelCore.UI.WriteLine("  MaxParallelThreads = " + ParallelThreads.ToString() + ".");
                }
                else
                    throw new System.Exception("Parallel must be 'true', 'false' or an integer >= 1.");
            }
            this.ThreadCount = ParallelThreads;
            FTimeStep = 1.0F / Timestep;
            if (! Names.TryGetParameter(Names.ClimateConfigFile, out var climateLibraryFileName))
            {
                ModelCore.UI.WriteLine($"  No ClimateConfigFile provided. Using climate files in ecoregion parameters: {Names.parameters["EcoregionParameters"].Value}.");
                ObservedClimate.Initialize();
            }
            PnETSpecies = new PnETSpecies();
            SpeciesParameters.LoadParameters(PnETSpecies);
            Hydrology.Initialize();
            SiteCohorts.Initialize();
            string PARunits = ((Parameter<string>)Names.GetParameter(Names.PARunits)).Value;
            if (PARunits != "umol" && PARunits != "W/m2")
                throw new Exception("PARunits are not 'umol' or 'W/m2'.");
            InitializeClimateLibrary(StartDate.Year); // John McNabb: initialize climate library after EcoregionPnET has been initialized
            Reproduction.SufficientResources = SufficientResources;
            Reproduction.Establish = Establish;
            Reproduction.AddNewCohort = AddNewCohort;
            Reproduction.MaturePresent = MaturePresent;
            Reproduction.PlantingEstablish = PlantingEstablish;
            SeedingAlgorithms SeedAlgorithm = (SeedingAlgorithms)Enum.Parse(typeof(SeedingAlgorithms), Names.parameters["SeedingAlgorithm"].Value);
            base.Initialize(ModelCore, SeedAlgorithm);
            ModelCore.UI.WriteLine("Spinning up biomass or reading from maps...");
            string InitialCommunitiesTxtFile = Names.GetParameter(Names.InitialCommunities).Value;
            string InitialCommunitiesMapFile = Names.GetParameter(Names.InitialCommunitiesMap).Value;
            InitialCommunitiesSpinup = Names.GetParameter(Names.InitialCommunitiesSpinup).Value;
            MinFolRatioFactor = ((Parameter<float>)Names.GetParameter(Names.MinFolRatioFactor, 0, float.MaxValue)).Value;
            Parameter<string> LitterMapFile;
            bool litterMapFile = Names.TryGetParameter(Names.LitterMap, out LitterMapFile);
            Parameter<string> WoodyDebrisMapFile;
            bool woodyDebrisMapFile = Names.TryGetParameter(Names.WoodyDebrisMap, out WoodyDebrisMapFile);
            InitializeSites(InitialCommunitiesTxtFile, InitialCommunitiesMapFile, ModelCore);
            if (litterMapFile)
                MapReader.ReadLitterFromMap(LitterMapFile.Value);
            if (woodyDebrisMapFile)
                MapReader.ReadWoodyDebrisFromMap(WoodyDebrisMapFile.Value);
            // Convert PnET cohorts to biomasscohorts
            foreach (ActiveSite site in ModelCore.Landscape)
            {
                SiteVars.UniversalCohorts[site] = SiteVars.SiteCohorts[site];
                if (SiteVars.SiteCohorts[site] != null && SiteVars.UniversalCohorts[site] == null)
                    throw new Exception("Cannot convert PnET SiteCohorts to biomass site cohorts");
            }
            ModelCore.RegisterSiteVar(SiteVars.UniversalCohorts, "Succession.UniversalCohorts");
            ISiteVar<SiteCohorts> PnETCohorts = ModelCore.Landscape.NewSiteVar<SiteCohorts>();
            foreach (ActiveSite site in ModelCore.Landscape)
            {
                PnETCohorts[site] = SiteVars.SiteCohorts[site];
                SiteVars.FineFuels[site] = SiteVars.Litter[site].Mass;
                IPnETEcoregionData ecoregion = PnETEcoregionData.GetPnETEcoregion(ModelCore.Ecoregion[site]);
                IHydrology hydrology = new Hydrology(ecoregion.FieldCap);
                float currentPressureHead = hydrology.PressureHeadTable.CalcSoilWaterPressureHead(hydrology.Water, ecoregion.SoilType);
                SiteVars.PressureHead[site] = currentPressureHead;
                SiteVars.FieldCapacity[site] = ecoregion.FieldCap / 10.0F; // cm volume (accounts for rooting depth)
                if (UsingClimateLibrary)
                {
                    SiteVars.ExtremeMinTemp[site] = (float)Enumerable
                        .Min(Climate.FutureEcoregionYearClimate[ecoregion.Index]
                        .Min(x => x.MonthlyTemp)) - (float)(3.0 * ecoregion.WinterSTD);
                    if (((Parameter<bool>)Names.GetParameter(Names.SoilIceDepth)).Value)
                    {
                        if (SiteVars.MonthlySoilTemp[site].Count() == 0)
                        {
                            // Calculations for soil temperature, ignoring snow?
                            // MG 20251009 -- replace parts with formulations in 
                            // PnET-Cohort Library, esp. soil thermal conductivity
                            float ThermalConductivity_theta = Soils.CalcThermalConductivitySoil_Watts(hydrology.SoilWaterContent, ecoregion.Porosity, ecoregion.SoilType) / Constants.Convert_kJperday_to_Watts;
                            float D = ThermalConductivity_theta / Hydrology_SaxtonRawls.GetCTheta(ecoregion.SoilType);  //m2/day
                            float Dmms = D * 1000000F / Constants.SecondsPerDay; // mm2/s
                            float d = (float)Math.Sqrt(2 * Dmms / Constants.omega);
                            float maxDepth = ecoregion.RootingDepth + ecoregion.LeakageFrostDepth;
                            float bottomFreezeDepth = maxDepth / 1000;
                            foreach (var year in Climate.SpinupEcoregionYearClimate[ecoregion.Index])
                            {
                                List<double> monthlyAirT = Climate.SpinupEcoregionYearClimate[ecoregion.Index][year.CalendarYear].MonthlyTemp;
                                double annualAirTemp = Climate.SpinupEcoregionYearClimate[ecoregion.Index][year.CalendarYear].MeanAnnualTemperature;
                                List<double> monthlyPrecip = Climate.SpinupEcoregionYearClimate[ecoregion.Index][year.CalendarYear].MonthlyPrecip;
                                SortedList<float, float> depthTempDict = new SortedList<float, float>();
                                SiteVars.MonthlyPressureHead[site] = new float[monthlyAirT.Count()];
                                SiteVars.MonthlySoilTemp[site] = new SortedList<float, float>[monthlyAirT.Count()];
                                for (int m = 0; m < monthlyAirT.Count(); m++)
                                {
                                    SiteVars.MonthlyPressureHead[site][m] = currentPressureHead;
                                    float DRz_snow = 1F; // Assume no snow in initialization
                                    // Damping ratio for moss - adapted from Kang et al. (2000) and Liang et al. (2014)
                                    float thermalDamping_Moss = (float)Math.Sqrt(2.0F * Constants.ThermalDiffusivityMoss / Constants.omega);
                                    float DRz_moss = (float)Math.Exp(-1.0F * ecoregion.MossDepth * thermalDamping_Moss);
                                    // Fill the tempDict with values
                                    float testDepth = 0F;
                                    float zTemp = 0F;
                                    int month = m + 1;
                                    int maxMonth = 0;
                                    int minMonth = 0;
                                    int mCount = 0;
                                    float tSum = 0F;
                                    float pSum = 0F;
                                    float tMax = float.MinValue;
                                    float tMin = float.MaxValue;
                                    if (m < 12)
                                    {
                                        mCount = Math.Min(12, monthlyAirT.Count());
                                        foreach (int z in Enumerable.Range(0, mCount))
                                        {
                                            tSum += (float)monthlyAirT[z];
                                            pSum += (float)monthlyPrecip[z];
                                            if (monthlyAirT[z] > tMax)
                                            {
                                                tMax = (float)monthlyAirT[z];
                                                maxMonth = z + 1;
                                            }
                                            if (monthlyAirT[z] < tMin)
                                            {
                                                tMin = (float)monthlyAirT[z];
                                                minMonth = z + 1;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        mCount = 12;
                                        foreach (int z in Enumerable.Range(m - 11, 12))
                                        {
                                            tSum += (float)monthlyAirT[z];
                                            pSum += (float)monthlyPrecip[z];
                                            if ((float)monthlyAirT[z] > tMax)
                                            {
                                                tMax = (float)monthlyAirT[z];
                                                maxMonth = month + z;
                                            }
                                            if ((float)monthlyAirT[z] < tMin)
                                            {
                                                tMin = (float)monthlyAirT[z];
                                                minMonth = month + z;
                                            }
                                        }
                                    }
                                    float annualTavg = tSum / mCount;
                                    float annualPcpAvg = pSum / mCount;
                                    float tAmplitude = (tMax - tMin) / 2;
                                    // Calculate depth to bottom of ice lens with FrostDepth
                                    while (testDepth <= bottomFreezeDepth)
                                    {
                                        float DRz = (float)Math.Exp(-1.0F * testDepth * d * ecoregion.FrostFactor); // adapted from Kang et al. (2000) and Liang et al. (2014); added FrostFactor for calibration
                                        // Calculate lag months from both max and min temperature months
                                        int lagMax = month + (3 - maxMonth);
                                        int lagMin = month + (minMonth - 5);
                                        if (minMonth >= 9)
                                            lagMin = month + (minMonth - 12 - 5);
                                        float lagAvg = (float)(lagMax + lagMin) / 2f;
                                        zTemp = (float)(annualAirTemp + tAmplitude * DRz_snow * DRz_moss * DRz * Math.Sin(Constants.omega * lagAvg - testDepth / d));
                                        depthTempDict[testDepth] = zTemp;
                                        if (testDepth == 0f)
                                            testDepth = 0.10f;
                                        else if (testDepth == 0.10f)
                                            testDepth = 0.25f;
                                        else
                                            testDepth += 0.25F;
                                    }
                                    SiteVars.MonthlySoilTemp[site][m] = Soils.CalculateMonthlySoilTemps(depthTempDict, ecoregion, 0, 0, hydrology, (float)monthlyAirT[m]);
                                }
                            }
                        }
                    }
                }
                else
                    SiteVars.ExtremeMinTemp[site] = 999;
            }
            ModelCore.RegisterSiteVar(PnETCohorts, "Succession.CohortsPnET");
        }

        /// <summary>
        /// This must be called *after* PnETPnETEcoregionData.Initialize() has been called
        /// </summary>
        private void InitializeClimateLibrary(int startYear = 0)
        {
            // John McNabb: initialize ClimateRegionData after initializing EcoregionPnet
            Parameter<string> climateLibraryFileName;
            UsingClimateLibrary = Names.TryGetParameter(Names.ClimateConfigFile, out climateLibraryFileName);
            if (UsingClimateLibrary)
            {
                ModelCore.UI.WriteLine($"Using climate library: {climateLibraryFileName.Value}.");
                Climate.Initialize(climateLibraryFileName.Value, false, ModelCore);
                ClimateRegionData.Initialize();
            }
            string PARunits = ((Parameter<string>)Names.GetParameter(Names.PARunits)).Value;
            if (PARunits == "umol")
                ModelCore.UI.WriteLine("Using PAR units of umol/m2/s.");
            else if (PARunits == "W/m2")
                ModelCore.UI.WriteLine("Using PAR units of W/m2.");
            else
                throw new ApplicationException(string.Format("PARunits units are not 'umol' or 'W/m2'"));
        }

        public void AddNewCohort(ISpecies species, ActiveSite site, string reproductionType, double biomassFrac = 1.0)
        {
            IPnETSpecies spc = PnETSpecies[species];
            bool addCohort = true;
            if (SiteVars.SiteCohorts[site].cohorts.ContainsKey(species))
            {
                // This should deliver only one KeyValuePair
                KeyValuePair<ISpecies, List<Cohort>> i = new List<KeyValuePair<ISpecies, List<Cohort>>>(SiteVars.SiteCohorts[site].cohorts.Where(o => o.Key == species))[0];
                List<Cohort> Cohorts = new List<Cohort>(i.Value.Where(o => o.Age < CohortBinSize));
                if (Cohorts.Count() > 0)
                    addCohort = false;
            }
            bool addSiteOutput = false;
            addSiteOutput = SiteOutputNames.ContainsKey(site) && addCohort;
            Cohort cohort = new Cohort(species, spc, (ushort)Date.Year, addSiteOutput ? SiteOutputNames[site] : null, biomassFrac, false);
            if (((Parameter<bool>)Names.GetParameter(Names.CohortStacking)).Value)
            {
                cohort.CanopyGrowingSpace = 1.0f;
                cohort.CanopyLayerProp = 1.0f;
            }
            addCohort = SiteVars.SiteCohorts[site].AddNewCohort(cohort);
            if (addCohort)
            {
                if (reproductionType == "plant")
                {
                    if (!SiteVars.SiteCohorts[site].SpeciesEstablishedByPlant.Contains(species))
                        SiteVars.SiteCohorts[site].SpeciesEstablishedByPlant.Add(species);
                }
                else if (reproductionType == "serotiny")
                {
                    if (!SiteVars.SiteCohorts[site].SpeciesEstablishedBySerotiny.Contains(species))
                        SiteVars.SiteCohorts[site].SpeciesEstablishedBySerotiny.Add(species);
                }
                else if (reproductionType == "resprout")
                {
                    if (!SiteVars.SiteCohorts[site].SpeciesEstablishedByResprout.Contains(species))
                        SiteVars.SiteCohorts[site].SpeciesEstablishedByResprout.Add(species);
                }
                else if (reproductionType == "seed")
                {
                    if (!SiteVars.SiteCohorts[site].SpeciesEstablishedBySeed.Contains(species))
                        SiteVars.SiteCohorts[site].SpeciesEstablishedBySeed.Add(species);
                }
            }
        }

        public bool MaturePresent(ISpecies species, ActiveSite site)
        {
            bool IsMaturePresent = SiteVars.SiteCohorts[site].IsMaturePresent(species);
            return IsMaturePresent;
        }

        protected override void InitializeSite(ActiveSite site)
        {
            lock (threadLock)
            {
                if (m == null)
                    m = new MyClock(ModelCore.Landscape.ActiveSiteCount);
                m.Next();
                m.WriteUpdate();
            }
            uint key = 0;
            allKeys.TryGetValue(site, out key);
            ICommunity initialCommunity = null;
            if (!sitesAndCommunities.TryGetValue(site, out initialCommunity))
                throw new ApplicationException(string.Format("Unable to retrieve initialCommunity for site: {0}", site.Location.Row + "," + site.Location.Column));
            if (!SiteCohorts.InitialSitesContainsKey(key))
            {
                // Create new sitecohorts from scratch
                SiteVars.SiteCohorts[site] = new SiteCohorts(StartDate, site, initialCommunity, UsingClimateLibrary, PlugIn.InitialCommunitiesSpinup, MinFolRatioFactor, SiteOutputNames.ContainsKey(site) ? SiteOutputNames[site] : null);
            }
            else
            {
                // Create new sitecohorts using initialcommunities data
                SiteVars.SiteCohorts[site] = new SiteCohorts(StartDate, site, initialCommunity, SiteOutputNames.ContainsKey(site) ? SiteOutputNames[site] : null);
            }
        }

        public override void InitializeSites(string initialCommunitiesText, string initialCommunitiesMap, ICore modelCore)
        {
            ModelCore.UI.WriteLine("   Loading initial communities from file \"{0}\" ...", initialCommunitiesText);
            DatasetParser parser = new DatasetParser(Timestep, modelCore.Species, additionalCohortParameters, initialCommunitiesText);
            IDataset communities = Landis.Data.Load<IDataset>(initialCommunitiesText, parser);
            List<ActiveSite> processFirst = new List<ActiveSite>();
            List<ActiveSite> processSecond = new List<ActiveSite>();
            ModelCore.UI.WriteLine("   Reading initial communities map \"{0}\" ...", initialCommunitiesMap);
            ProcessInitialCommunitiesMap(initialCommunitiesMap, communities, ref processFirst, ref processSecond);
            if (this.ThreadCount != 1)
            {
                // Handle creation of initial community sites first
                Parallel.ForEach
                (
                    processFirst, new ParallelOptions
                    {
                        MaxDegreeOfParallelism = this.ThreadCount
                    }, site =>
                    {
                        InitializeSite(site);
                    }
                );
                Parallel.ForEach
                (
                    processSecond, new ParallelOptions
                    {
                        MaxDegreeOfParallelism = this.ThreadCount
                    }, site =>
                    {
                        InitializeSite(site);
                    }
                );
            }
            else
            {
                // First, process sites so that the initial communities are set up
                foreach (ActiveSite site in processFirst)
                    InitializeSite(site);
                foreach (ActiveSite site in processSecond)
                    InitializeSite((ActiveSite)site);
            }
        }

        protected override void AgeCohorts(ActiveSite site,
                                           ushort years,
                                           int? successionTimestep)
        {
            // Date starts at 1/15/Year
            DateTime date = new DateTime(StartDate.Year + ModelCore.CurrentTime - Timestep, 1, 15);
            DateTime EndDate = date.AddYears(years);
            IPnETEcoregionData PnETEcoregion = PnETEcoregionData.GetPnETEcoregion(ModelCore.Ecoregion[site]);
            List<IPnETEcoregionVars> climate_vars = UsingClimateLibrary ? PnETEcoregionData.GetClimateRegionData(PnETEcoregion, date, EndDate) : PnETEcoregionData.GetData(PnETEcoregion, date, EndDate);
            SiteVars.SiteCohorts[site].Grow(climate_vars);
            SiteVars.SiteCohorts[site].DisturbanceTypesReduced.Clear();
            Date = EndDate;
        }

        // Required function - not used within PnET-Succession
        public override byte ComputeShade(ActiveSite site)
        {
            return 0;
        }

        public override void Run()
        {
            if (Timestep > 0)
                ClimateRegionData.SetAllEcoregionsFutureAnnualClimate(ModelCore.CurrentTime);
            base.Run();
        }

        // This is a Delegate method to base succession.
        // Not used within PnET-Succession
        public bool SufficientResources(ISpecies species, ActiveSite site)
        {
            return true;
        }

        /// <summary>
        /// Determines if a species can establish on a site.
        /// This is a Delegate method to base succession.
        /// </summary>
        public bool Establish(ISpecies species, ActiveSite site)
        {
            IPnETSpecies spc = PnETSpecies[species];
            bool Establish = SiteVars.SiteCohorts[site].EstablishmentProbability.HasEstablished(spc);
            return Establish;
        }

        /// <summary>
        /// Determines if a species can be planted on a site 
        /// (all conditions are satisfied).
        /// This is a Delegate method to base succession.
        /// </summary>
        public bool PlantingEstablish(ISpecies species, ActiveSite site)
        {
            return true;
        }

        /// <summary>
        /// Reads the initial communities map, finds all unique site keys, and sets aside sites to process first and second
        /// </summary>
        private void ProcessInitialCommunitiesMap(string initialCommunitiesMap, 
                                                  IDataset communities,
                                                  ref List<ActiveSite> processFirst,
                                                  ref List<ActiveSite> processSecond)
        {
            IInputRaster<UIntPixel> map = ModelCore.OpenRaster<UIntPixel>(initialCommunitiesMap);
            Dictionary<uint, ActiveSite> uniqueKeys = new Dictionary<uint, ActiveSite>();
            using (map)
            {
                UIntPixel pixel = map.BufferPixel;
                foreach (Site site in ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    uint mapCode = pixel.MapCode.Value;
                    if (! site.IsActive)
                        continue;
                    ActiveSite activeSite = (ActiveSite)site;
                    var initialCommunity = communities.Find(mapCode);
                    if (initialCommunity == null)
                        throw new ApplicationException(string.Format("Unknown map code for initial community: {0}", mapCode));
                    sitesAndCommunities.Add(activeSite, initialCommunity);
                    uint key = SiteCohorts.ComputeKey((ushort)initialCommunity.MapCode, Globals.ModelCore.Ecoregion[site].MapCode);
                    if (! uniqueKeys.ContainsKey(key))
                    {
                        uniqueKeys.Add(key, activeSite);
                        processFirst.Add(activeSite);
                    }
                    else
                        processSecond.Add(activeSite);
                    if (!allKeys.ContainsKey(activeSite))
                        allKeys.Add(activeSite, key);
                }
            }
        }

        public override void AddCohortData()
        {
            // CUSTOM DYNAMIC PARAMETERS GO HERE
            return;
        }
    }
}
