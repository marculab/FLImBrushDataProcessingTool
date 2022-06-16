ini = IniConfig();
ini.AddSections({'Software Version'});
ini.AddKeys('Software Version', 'MAJOR', 1);
ini.AddKeys('Software Version', 'MINOR', 0);
ini.AddKeys('Software Version', 'PATCH', 0);
[filepath,name,ext] = fileparts(pwd);
ini.WriteFile(fullfile(filepath,'SoftwareVersion.ini'));