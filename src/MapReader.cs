// Original Authors: Austen Ruzicka and Robert Scheller

// NOTE: DoublePixel --> Landis.SpatialModeling
// NOTE: IInputRaster --> Landis.SpatialModeling
// NOTE: InputValueException --> Landis.Utilities
// NOTE: Site --> Landis.SpatialModeling
// NOTE: SiteVars --> Landis.Library.PnETCohorts

using System;
using System.IO;
// using Landis.Core;
using Landis.Library.PnETCohorts;
using Landis.SpatialModeling;using Landis.Utilities;

namespace Landis.Extension.Succession.BiomassPnET
{
    /// <summary>
    /// Methods to read maps instead of using domain spin-up
    /// </summary>
    public static class MapReader
    {
        static private double maxLeafLitter = 5176.124;
        static private double maxWoodDebris = 95007.64;
        static private double minLeafLitter = 0;
        static private double minWoodDebris = 0;

        public static void ReadWoodDebrisFromMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path);
            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    int mapValue = (int)pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < minWoodDebris || mapValue > maxWoodDebris)
                            throw new InputValueException(mapValue.ToString(),
                                                          "Wood debris value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, minWoodDebris, maxWoodDebris, site.Location.Row, site.Location.Column);
                        SiteVars.WoodDebris[site].InitialMass = mapValue;
                        SiteVars.WoodDebris[site].Mass = mapValue;
                    }
                }
            }
        }

        public static void ReadLeafLitterFromMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path);
            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    int mapValue = (int)pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < minLeafLitter || mapValue > maxLeafLitter)
                            throw new InputValueException(mapValue.ToString(),
                                                          "Leaf litter value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, minLeafLitter, maxLeafLitter, site.Location.Row, site.Location.Column);
                        SiteVars.LeafLitter[site].InitialMass = mapValue;
                        SiteVars.LeafLitter[site].Mass = mapValue;
                    }
                }
            }
        }

        private static IInputRaster<DoublePixel> MakeDoubleMap(string path)
        {
            PlugIn.ModelCore.UI.WriteLine("  Read in data from {0}", path);
            IInputRaster<DoublePixel> map;
            try
            {
                map = PlugIn.ModelCore.OpenRaster<DoublePixel>(path);
            }
            catch (FileNotFoundException)
            {
                string msg = string.Format("Error: The file {0} does not exist", path);
                throw new ApplicationException(msg);
            }
            if (map.Dimensions != PlugIn.ModelCore.Landscape.Dimensions)
            {
                string msg = string.Format("Error: The input map {0} does not have the same dimension (row, column) as the scenario ecoregions map", path);
                throw new ApplicationException(msg);
            }
            return map;
        }
    }
}
