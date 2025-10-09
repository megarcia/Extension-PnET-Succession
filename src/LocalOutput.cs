using System.Collections.Generic;
using System.IO;

namespace Landis.Extension.Succession.BiomassPnET 
{
    public class LocalOutput
    {
        public static string PnETOutputSites;
        private List<string> FileContent;
        public string FileName { get; private set; }
        public string SiteName { get; private set; }
        public string Path { get; private set; }

        public LocalOutput(string SiteName, string FileName, string Header)
        {
            this.SiteName = SiteName;
            Path = "Output" + Path.DirectorySeparatorChar + PnETOutputSites + Path.DirectorySeparatorChar + SiteName + Path.DirectorySeparatorChar;
            this.FileName = FileName;
            if (File.Exists(Path + FileName))
                File.Delete(Path + FileName);
            if (Directory.Exists(Path) == false)
                Directory.CreateDirectory(Path);
            FileContent = new List<string>(new string[] { Header });
            Write();
        }

        public void Add(string s)
        {
            FileContent.Add(s);
        }

        public void Write()
        {
            try
            {
                StreamWriter sw = new StreamWriter(Path.Combine(Path, FileName), true);
                foreach (string line in FileContent)
                    sw.WriteLine(line);
                sw.Close();
                FileContent.Clear();
                return;
            }
            catch (IOException e)
            {
                PlugIn.ModelCore.UI.WriteLine("Cannot write to " + Path.Combine(Path, FileName) + " " + e.Message);
            }
        }
    }
}
