function E = w2f_machine_info()
%W2F_MACHINE_INFO Capture reproducibility metadata without requiring toolboxes.
E=struct();
E.MATLABVersion=version;
E.Computer=computer;
E.Hostname=getenv('COMPUTERNAME');
if isempty(E.Hostname),E.Hostname=getenv('HOSTNAME');end
E.ProcessorIdentifier=getenv('PROCESSOR_IDENTIFIER');
E.NumberOfProcessors=getenv('NUMBER_OF_PROCESSORS');
E.OS=getenv('OS');
try,E.NumCoresFeature=feature('numcores');catch,E.NumCoresFeature=NaN;end
try
    [st,txt]=system('wmic cpu get Name,NumberOfCores,NumberOfLogicalProcessors /value');
    if st~=0 || isempty(strtrim(txt))
        [st,txt]=system('powershell -NoProfile -Command "Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors | Format-List"');
    end
    if st==0,E.CPUQuery=strtrim(txt);else,E.CPUQuery='unavailable';end
catch
    E.CPUQuery='unavailable';
end
try
    E.JavaVersion=version('-java');
catch
    E.JavaVersion='unavailable';
end
end
