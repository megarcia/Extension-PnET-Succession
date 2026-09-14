using System.Collections.Generic;
using System.IO;

namespace Landis.Extension.Succession.BiomassPnET 
{
    public class LocalOutput
    {
        public static string PNEToutputsites;
        private List<string> FileContent;
        public string FileName { get; private set; }
        public string SiteName { get; private set; }
        public string FilePath { get; private set; }

        public LocalOutput(string SiteName, string FileName, string Header)
        {
            this.SiteName = SiteName;
            this.FileName = FileName;
            FilePath = "Output" + Path.DirectorySeparatorChar + PNEToutputsites + Path.DirectorySeparatorChar + SiteName + Path.DirectorySeparatorChar;
            if (File.Exists(FilePath + FileName))
                File.Delete(FilePath + FileName);
            if (Directory.Exists(FilePath) == false)
                Directory.CreateDirectory(FilePath);
            FileContent = new List<string>(new string[] { Header });
            Write();
        }

        public void Add(string s)
        {
            FileContent.Add(s);
        }

        public void Write()
        {
            while (true)
            {
                try
                {
                    StreamWriter sw = new StreamWriter(Path.Combine(FilePath, FileName), true);
                    foreach (string line in FileContent)
                        sw.WriteLine(line);
                    sw.Close();
                    FileContent.Clear();
                    return;
                }
                catch (IOException e)
                {
                    PlugIn.ModelCore.UI.WriteLine("Cannot write to " + Path.Combine(FilePath, FileName) + " " + e.Message);
                }
            }
        }
    }
}
