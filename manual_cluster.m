function varargout = manual_cluster(varargin)
% MANUAL_CLUSTER MATLAB code for manual_cluster.fig
%      MANUAL_CLUSTER, by itself, creates a new MANUAL_CLUSTER or raises the existing
%      singleton*.
%
%      H = MANUAL_CLUSTER returns the handle to a new MANUAL_CLUSTER or the handle to
%      the existing singleton*.
%
%      MANUAL_CLUSTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUAL_CLUSTER.M with the given input arguments.
%
%      MANUAL_CLUSTER('Property','Value',...) creates a new MANUAL_CLUSTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manual_cluster_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manual_cluster_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manual_cluster

% Last Modified by GUIDE v2.5 03-Feb-2015 19:55:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manual_cluster_OpeningFcn, ...
                   'gui_OutputFcn',  @manual_cluster_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before manual_cluster is made visible.
function manual_cluster_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manual_cluster (see VARARGIN)

% Choose default command line output for manual_cluster
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes manual_cluster wait for user response (see UIRESUME)
% uiwait(handles.figure1);

spikemat = varargin{1};
spkid = varargin{2};
spktime = varargin{3};
spcorr=spikemat(:,:)*spikemat(:,:)';
[u,~,~]=svd(spcorr);
nc=2; %only for two dimensional
proj_u=spikemat(:,:)'*u(:,1:nc);
axes(handles.axes_cluster_plot)
plot(proj_u(:,1),proj_u(:,2),'.r')
idx=zeros(size(spikemat,2),1);
clusters=[]; spcounts=[]; centers=[]; spmat=[];
handles.spcounts=spcounts;
guidata(hObject, handles);
handles.centers=centers;
guidata(hObject, handles);
handles.spmat=spmat;
guidata(hObject, handles);
save('spikedata.mat','proj_u','idx','clusters','u','spcounts','spmat','centers','spikemat','spkid','spktime','-append');

uiwait(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = manual_cluster_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.spcounts;  
varargout{2} = handles.centers;
varargout{3} = handles.spmat;
 delete(hObject);
 
% --- Executes on button press in pushbutton_addcluster.
function pushbutton_addcluster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('spikedata.mat','proj_u','idx','clusters');
popup_label=[];
cell_selectcluster = get(handles.popupmenu_selectcluster,'string');
length_cell_selectcluster = length(cell_selectcluster);
eval(sprintf('popup_label=''cluster%i'';',max(idx)+1))
cell_selectcluster{length_cell_selectcluster + 1} = popup_label;
set(handles.popupmenu_selectcluster, 'String', cell_selectcluster);
set(handles.popupmenu_selectcluster, 'value', length(cell_selectcluster));



convex_boundary_x=[]; convex_boundary_y=[]; k=[];
[convex_boundary_x,convex_boundary_y] = ginput();
convex_boundary_x = [convex_boundary_x;convex_boundary_x(1)];
convex_boundary_y = [convex_boundary_y;convex_boundary_y(1)];
k = convhull(convex_boundary_x,convex_boundary_y);
convex_boundary_x = convex_boundary_x(k);
convex_boundary_y = convex_boundary_y(k);
eval(sprintf('clusters.cluster%i.convex_boundary(:,1)=convex_boundary_x;',max(idx)+1))
eval(sprintf('clusters.cluster%i.convex_boundary(:,2)=convex_boundary_y;',max(idx)+1))
in = inpolygon(proj_u(:,1),proj_u(:,2),convex_boundary_x,convex_boundary_y);
idx(in)=max(idx)+1;

cluster_color=hsv(10);
axes(handles.axes_cluster_plot)
for i=1:max(idx)+1
    t=[];
    eval(sprintf('t=proj_u(find(idx==%i),:);',i-1))
    if ~isempty(t)
    plot(t(:,1),t(:,2),'.','color',cluster_color(i,:))
    end
    hold on
end
save('spikedata.mat','idx','clusters','-append');
hold off




% --- Executes on button press in pushbutton_removecluster.
function pushbutton_removecluster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removecluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load('spikedata.mat','proj_u','idx','clusters');

cell_selectcluster = cellstr(get(handles.popupmenu_selectcluster,'string'));
selected_cluster=get(handles.popupmenu_selectcluster, 'value');
popup_label=char(cell_selectcluster(selected_cluster));
cluster_no=str2num(popup_label( regexp(popup_label,'\d')));
eval(sprintf('clusters.%s=[];',popup_label))
cell_selectcluster(selected_cluster) =[];
set(handles.popupmenu_selectcluster, 'String', cell_selectcluster);
if selected_cluster == 1
    set(handles.popupmenu_selectcluster, 'value', length(cell_selectcluster));
else
set(handles.popupmenu_selectcluster, 'value', selected_cluster-1);
end

ax= idx==cluster_no;
idx(ax)=0;

cluster_color=hsv(10);
axes(handles.axes_cluster_plot)
for i=1:max(idx)+1
    eval(sprintf('t=proj_u(find(idx==%i),:);',i-1))
    plot(t(:,1),t(:,2),'.','color',cluster_color(i,:))
    hold on
end
save('spikedata.mat','idx','clusters','-append');
hold off


% --- Executes on button press in pushbutton_recluster.
function pushbutton_recluster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_recluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load('spikedata.mat','proj_u','idx','clusters');
cell_selectcluster = get(handles.popupmenu_selectcluster,'string');
selected_cluster=get(handles.popupmenu_selectcluster, 'value');
popup_label=char(cell_selectcluster(selected_cluster));
eval(sprintf('clusters.%s=[];',popup_label))
cluster_no=str2num(popup_label( regexp(popup_label,'\d')));
ax= idx==cluster_no; % getting index of cluster to be reclustered
idx(ax)=0; % resetting the selected cluster to belong to overall raw cluster

% RECLUSTERING
convex_boundary_x=[]; convex_boundary_y=[]; k=[];

[convex_boundary_x,convex_boundary_y] = ginput();
convex_boundary_x = [convex_boundary_x;convex_boundary_x(1)]; 
convex_boundary_y = [convex_boundary_y;convex_boundary_y(1)];
	k = convhull(convex_boundary_x,convex_boundary_y); 
	convex_boundary_x = convex_boundary_x(k); 
	convex_boundary_y = convex_boundary_y(k);
in = inpolygon(proj_u(:,1),proj_u(:,2),convex_boundary_x,convex_boundary_y);
idx(in)=cluster_no; % alloted the same cluster number

eval(sprintf('clusters.cluster%i.convex_boundary(:,1)=convex_boundary_x;',cluster_no))
eval(sprintf('clusters.cluster%i.convex_boundary(:,2)=convex_boundary_y;',cluster_no))
cluster_color=hsv(10);
axes(handles.axes_cluster_plot)
for i=1:max(idx)+1
    t=[];
    eval(sprintf('t=proj_u(find(idx==%i),:);',i-1))
    if ~isempty(t)
    plot(t(:,1),t(:,2),'.','color',cluster_color(i,:))
    end
    hold on
end
save('spikedata.mat','idx','clusters','-append');
hold off



% --- Executes on selection change in popupmenu_selectcluster.
function popupmenu_selectcluster_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_selectcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_selectcluster contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_selectcluster

load('spikedata.mat','proj_u','idx','clusters');
contents = cellstr(get(hObject,'String'));
selected_cluster=char(contents(get(hObject,'value')));
eval(sprintf('boundary=clusters.%s.convex_boundary;',selected_cluster))
cluster_color=hsv(10);
axes(handles.axes_cluster_plot)
for i=1:max(idx)+1
    t=[];
    eval(sprintf('t=proj_u(find(idx==%i),:);',i-1))
    if ~isempty(t)
    plot(t(:,1),t(:,2),'.','color',cluster_color(i,:))
    end
    hold on
end
hold on
plot(boundary(:,1),boundary(:,2))
hold off


% --- Executes during object creation, after setting all properties.
function popupmenu_selectcluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_selectcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton_loadcluster.
function togglebutton_loadcluster_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_loadcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_loadcluster

load('spikedata.mat','spktime','spkid','idx','u','proj_u','spikemat','spcounts','centers','spmat');
cell_selectcluster = get(handles.popupmenu_selectcluster,'string');


for j=1:length(fields(spkid))
    eval(sprintf('itersids=spkid.iter%i(:);',j))
    for sn=1:length(itersids)
        
        for i=1:length(cell_selectcluster)
            ad=[];
            popup_label=char(cell_selectcluster(i));
            cluster_no=str2num(popup_label( regexp(popup_label,'\d')));
            ad=idx==cluster_no;
            spike_data=(spikemat(:,ad));
            if ad(itersids(sn))
                iterclass(sn)=i;
            end
        end
    end
    if exist('iterclass')
        eval(sprintf('spkiterclass.iter%i=iterclass;',j))
        clear iterclass
    else
        eval(sprintf('spkiterclass.iter%i=[];',j))
    end
    
end

for unit=1:1:length(cell_selectcluster)
    for kk=1:length(fields(spkid))
        eval(sprintf('spcounts.cl%i.iter%i=spktime.iter%i(find(spkiterclass.iter%i==unit));',unit,kk,kk,kk));
    end
end

         

for i=1:length(cell_selectcluster)
    ad=[]; proj=[];
    popup_label=char(cell_selectcluster(i));
    cluster_no=str2num(popup_label( regexp(popup_label,'\d')));
    ad=idx==cluster_no;
    proj=mean(proj_u((ad),:));
    centers(:,i)=u(:,1:2)*proj';
end

for ii=1:length(cell_selectcluster)
            ad=[]; spike_data=[];
            popup_label=char(cell_selectcluster(ii));
            cluster_no=str2num(popup_label( regexp(popup_label,'\d')));
            ad=idx==cluster_no;
            spike_data=(spikemat(:,ad));
            eval(sprintf('spmat.cl%i=spike_data;',ii));
end
% save('spikedata.mat','spcounts','centers','-append');
handles.spcounts=spcounts; 
guidata(hObject, handles); 
handles.centers=centers;
guidata(hObject, handles); 
handles.spmat=spmat; 
guidata(hObject, handles); 

 close(handles.figure1);
 
% --- Executes on button press in pushbutton_waveform.
function pushbutton_waveform_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_waveform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load('spikedata.mat','spikemat','idx','u','proj_u');
cell_selectcluster = get(handles.popupmenu_selectcluster,'string');
selected_cluster=get(handles.popupmenu_selectcluster, 'value');
popup_label=char(cell_selectcluster(selected_cluster));
cluster_no=str2num(popup_label( regexp(popup_label,'\d')));
ad= idx==cluster_no;
proj=mean(proj_u((ad),:));
centers(:,1)=u(:,1:2)*proj';
spike_data=(spikemat(:,ad));
find_figure('waveform');
plot(spike_data,'color',[0.7 0.7 0.7])
hold on
plot(centers,'r')
hold off








% --- Executes on button press in pushbutton_isi.
function pushbutton_isi_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_isi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_autocorr.
function pushbutton_autocorr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_autocorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure


if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, call UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

