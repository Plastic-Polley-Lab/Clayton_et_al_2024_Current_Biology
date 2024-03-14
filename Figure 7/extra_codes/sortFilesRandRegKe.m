function struct_tosca=sortFilesRandRegKe(runFile)
% SORTFILESRANDREG read the tosca data of a run of rand reg trials and
% create a structure registering the set and cycle size of each trial.

[Data Params]=tosca_read_run(runFile);

for iTrial=1:length(Data)
    data=Data{iTrial}.RandReg.Speaker.File;
    trial=tosca_read_trial(Params, Data, iTrial);
    
    %Finds the onset and offset of the state during which the stimulus is
    %played
    idxStateChange=find((trial.State_Change(2:end)-trial.State_Change(1:end-1))==1)+1;
    tOnset=trial.Time_s(idxStateChange(1))-trial.Time_s(1);
    tOffset=trial.Time_s(idxStateChange(2))-trial.Time_s(1);
    %Set #9 corresponds to the RAND-RAND control. In the output, control is
    %set#4 (indicated in the "type" subfield).
%     if data.scale==9
%         data.set=4;
%         struct_tosca.set(data.set).type='RAND-RAND';
%     else
%          struct_tosca.set(data.set).type='RAND-REG';
%     end
    
  
    %Rearrage by cycle size (cyc=[4 6 8 10 12] mapped to iCyc=[1 2 3 4 5])
%     iCyc=(data.cyc-2)/2;
    struct_tosca(iTrial).rand=data.rand;
    struct_tosca(iTrial).scale=data.scale;
    struct_tosca(iTrial).tOnset=tOnset;
    struct_tosca(iTrial).tOffset=tOffset;
    struct_tosca(iTrial).block=Data{iTrial}.block;
    struct_tosca(iTrial).trial=Data{iTrial}.trial;
    struct_tosca(iTrial).N=Data{iTrial}.N;

end
end