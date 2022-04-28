function varargout = spike_removal_gui(varargin)
% SPIKE_REMOVAL_GUI MATLAB code for spike_removal_gui.fig
%      SPIKE_REMOVAL_GUI, by itself, creates a new SPIKE_REMOVAL_GUI or raises the existing
%      singleton*.
%
%      H = SPIKE_REMOVAL_GUI returns the handle to a new SPIKE_REMOVAL_GUI or the handle to
%      the existing singleton*.
%
%      SPIKE_REMOVAL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPIKE_REMOVAL_GUI.M with the given input arguments.
%
%      SPIKE_REMOVAL_GUI('Property','Value',...) creates a new SPIKE_REMOVAL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spike_removal_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spike_removal_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spike_removal_gui

% Last Modified by GUIDE v2.5 22-Jul-2015 12:26:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spike_removal_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @spike_removal_gui_OutputFcn, ...
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


% --- Executes just before spike_removal_gui is made visible.
function spike_removal_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spike_removal_gui (see VARARGIN)

% Choose default command line output for spike_removal_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spike_removal_gui wait for user response (see UIRESUME)
% uiwait(handles.spike_removal_gui);

spikemat = varargin{1};
spkid=varargin{2};
spktime=varargin{3};
idx = varargin{4};
%we=varargin{5};



handles.spikemat=spikemat;
guidata(hObject, handles);
handles.spkid=spkid;
guidata(hObject, handles);
handles.spktime=spktime;
guidata(hObject, handles);
handles.idx=idx; 
guidata(hObject, handles);


spcorr=spikemat*spikemat';
[u,~,~]=svd(spcorr);
nc=3; %only for two dimensional
proj_u=spikemat(:,:)'*u(:,1:nc);


axes(handles.axes_cluster_plot)
axes_size(1)=min(proj_u(:,1))- 0.1 * min(proj_u(:,1));
axes_size(3)=max(proj_u(:,1))+ 0.1 * max(proj_u(:,1));
axes_size(2)=min(proj_u(:,2))- 0.1 * min(proj_u(:,2));
axes_size(4)=max(proj_u(:,2))+ 0.1 * max(proj_u(:,2));

% if we==1
%     dm=[1 2];
% elseif we==2
%     dm=[1 3];
% elseif we==3
%     dm=[2 3];
% end
dm=[1 2];
in_empty=idx==-1;

axes(handles.axes_cluster_plot)
        
        if any(in_empty)
            plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
        else
            plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
        end
set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);

 axes(handles.axes_waveform_plot)
        
        if any(in_empty)
           plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
        else
           plot(spikemat,'color',[0.7 0.7 .7]) 
        end
set(gca, 'XLim', [0 size(spikemat,1)]);



clusters=[]; in=[];
proj_u_copy=proj_u;
spikemat_copy=spikemat;
save('spikedata_removal.mat','proj_u','idx','clusters','u','spikemat','spkid','spktime','in','spikemat_copy','proj_u_copy','axes_size','-append');

uiwait


% --- Outputs from this function are returned to the command line.
function varargout = spike_removal_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


varargout{1} = handles.spikemat;
varargout{2} = handles.spkid;
varargout{3} = handles.spktime;
varargout{4} = handles.idx;
delete(hObject); 


% uiresume(spike_removal_gui)


% --- Executes on selection change in popupmenu_selection_method.
function popupmenu_selection_method_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_selection_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_selection_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_selection_method

load('spikedata_removal.mat','proj_u','idx','axes_size','spikemat');

% if we==1
%     dm=[1 2];
% elseif we==2
%     dm=[1 3];
% elseif we==3
%     dm=[2 3];
% end
dm=[1 2];
selected_method=get(handles.popupmenu_selection_method, 'value');

  in=[];
 in_empty=idx==-1; % specifying the spike already removed by user
% if we==1
%     dm=[1 2];
% elseif we==2
%     dm=[1 3];
% elseif we==3
%     dm=[2 3];
% end
dm=[1 2];
switch selected_method
    
    case 1
        
       
        axes(handles.axes_cluster_plot)
axes_size(1)=min(proj_u(:,1))- 0.1 * min(proj_u(:,1));
axes_size(3)=max(proj_u(:,1))+ 0.1 * max(proj_u(:,1));
axes_size(2)=min(proj_u(:,2))- 0.1 * min(proj_u(:,2));
axes_size(4)=max(proj_u(:,2))+ 0.1 * max(proj_u(:,2));
        
        if any(in_empty)
            plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
        else
            plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
        end
        
        set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
        
        axes(handles.axes_waveform_plot)
        
        if any(in_empty)
           plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
        else
           plot(spikemat,'color',[0.7 0.7 .7]) 
        end
        
        set(gca, 'XLim', [0 size(spikemat,1)]);
        
    case 2
        
       axes(handles.axes_cluster_plot)
axes_size(1)=min(proj_u(:,1))- 0.1 * min(proj_u(:,1));
axes_size(3)=max(proj_u(:,1))+ 0.1 * max(proj_u(:,1));
axes_size(2)=min(proj_u(:,2))- 0.1 * min(proj_u(:,2));
axes_size(4)=max(proj_u(:,2))+ 0.1 * max(proj_u(:,2));
        
        if any(in_empty)
            plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
        else
            plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
        end
          set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);

        
        axes(handles.axes_waveform_plot)
        
        if any(in_empty)
           plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
        else
           plot(spikemat,'color',[0.7 0.7 .7]) 
        end
        set(gca, 'XLim', [0 size(spikemat,1)]);
        
        X=[]; Y=[]; selected_points=[];  spikes_idx=[];
        
        axes(handles.axes_cluster_plot)
        [X,Y]= getpts;
        selected_points=[X,Y];
        spikes_idx = knnsearch(proj_u,selected_points);
        
        in=zeros(1,size(spikemat,2));
        in(spikes_idx)=1;
        
        axes(handles.axes_cluster_plot)
        hold on
        plot(proj_u(spikes_idx,1),proj_u(spikes_idx,2),'*')
        hold off
         set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);

        axes(handles.axes_waveform_plot)
        hold on
        plot(spikemat(:,spikes_idx))
        hold off
        set(gca, 'XLim', [0 size(spikemat,1)]);

        save('spikedata_removal.mat','in','-append');
        
        set(handles.popupmenu_cluster_select, 'value', 1);  
       
    case 3
        
       axes(handles.axes_cluster_plot)
axes_size(1)=min(proj_u(:,1))- 0.1 * min(proj_u(:,1));
axes_size(3)=max(proj_u(:,1))+ 0.1 * max(proj_u(:,1));
axes_size(2)=min(proj_u(:,2))- 0.1 * min(proj_u(:,2));
axes_size(4)=max(proj_u(:,2))+ 0.1 * max(proj_u(:,2));
        if any(in_empty)
           plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
        else
           plot(spikemat,'color',[0.7 0.7 .7]) 
        end
        set(gca, 'XLim', [0 size(spikemat,1)]);

        axes(handles.axes_cluster_plot)
        
        if any(in_empty)
            plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
        else
            plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
        end
            set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
        
        X=[]; Y=[];
        axes(handles.axes_waveform_plot)
        [X,Y]= getpts;
%         in_mat=false(size(spikemat,2),length(X));
        
        for k=1:length(X)
            part_data=[];
            
          if floor(X(k))~=X(k) && floor(X(k))~=0
             part_data=spikemat(floor(X(k)):ceil(X(k)),:);
          elseif floor(X(k))== 0
             part_data=spikemat(1,:); 
          else
             part_data=spikemat(X(k),:); 
          end 
          
%             for i=1:size(part_data,2)
                
                spikes_idx(k) = knnsearch(part_data',repmat(Y(k),1,size(part_data,1)));
                
%                 if size(part_data,1)==2   
%                 in_mat(i,k)=any(part_data(1,i)>=Y(k)) & any(part_data(2,i)<=Y(k));
%                 else
%                     pr_spike=find(part_data==Y(k));
%                     in_mat(pr_spike,k)=1;
%                 end
%             end
        end
        
         in=zeros(1,size(spikemat,2));
         in(spikes_idx)=1;
        
        axes(handles.axes_waveform_plot)
        hold on
        plot(spikemat(:,spikes_idx))
        hold off
        set(gca, 'XLim', [0 size(spikemat,1)]);
        
        axes(handles.axes_cluster_plot)
        hold on
        plot(proj_u(spikes_idx,1),proj_u(spikes_idx,2),'*')
        hold off
             set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
        
        save('spikedata_removal.mat','in','-append');
        
        set(handles.popupmenu_cluster_select, 'value', 1); 
        
    case 4
        
        rect=[];
        set(handles.popupmenu_cluster_select, 'value', 1); 
        
        axes(handles.axes_waveform_plot)
        axes(handles.axes_cluster_plot)
axes_size(1)=min(proj_u(:,1))- 0.1 * min(proj_u(:,1));
axes_size(3)=max(proj_u(:,1))+ 0.1 * max(proj_u(:,1));
axes_size(2)=min(proj_u(:,2))- 0.1 * min(proj_u(:,2));
axes_size(4)=max(proj_u(:,2))+ 0.1 * max(proj_u(:,2));
        
        if any(in_empty)
           plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
        else
           plot(spikemat,'color',[0.7 0.7 .7]) 
        end
        set(gca, 'XLim', [0 size(spikemat,1)]);
        
        axes(handles.axes_cluster_plot)
        
        if any(in_empty)
            plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
        else
            plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
        end
           set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
        
        rect=getrect(handles.axes_waveform_plot);
        if rect(1)>=1
        part_data=spikemat(floor(rect(1)):ceil(rect(1)+rect(3)),:);
        else
        part_data=spikemat(1:ceil(rect(1)+rect(3)),:);   
        end
        for i=1:size(part_data,2)
            in(i)=any(part_data(:,i)>=rect(2)) & any(part_data(:,i)<=rect(2)+rect(4));
        end
        
        axes(handles.axes_waveform_plot)
        hold on
        plot(spikemat(:,logical(in)),'r')
        hold off
        set(gca, 'XLim', [0 size(spikemat,1)]);
        
        axes(handles.axes_cluster_plot)
        hold on
        plot(proj_u(logical(in),1),proj_u(logical(in),2),'.r')
        hold off
            set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
            
       save('spikedata_removal.mat','in','-append');
       
       
    case 5
        
        set(handles.popupmenu_cluster_select, 'value', 1); 
        
        axes(handles.axes_cluster_plot)
        
        if any(in_empty)
            plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
        else
            plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
        end
            set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
        
        axes(handles.axes_waveform_plot)
        
        if any(in_empty)
           plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
        else
           plot(spikemat,'color',[0.7 0.7 .7]) 
        end
        set(gca, 'XLim', [0 size(spikemat,1)]);

        axes(handles.axes_cluster_plot)
        
        convex_boundary_x=[]; convex_boundary_y=[]; k=[];
        [convex_boundary_x,convex_boundary_y] = ginput();
        convex_boundary_x = [convex_boundary_x;convex_boundary_x(1)];
        convex_boundary_y = [convex_boundary_y;convex_boundary_y(1)];
        k = convhull(convex_boundary_x,convex_boundary_y);
        convex_boundary_x = convex_boundary_x(k);
        convex_boundary_y = convex_boundary_y(k);
        in = inpolygon(proj_u(:,1),proj_u(:,2),convex_boundary_x,convex_boundary_y);
        
        axes(handles.axes_waveform_plot)
        hold on
        plot(spikemat(:,logical(in)),'r')
        hold off
        set(gca, 'XLim', [0 size(spikemat,1)]);
        
        axes(handles.axes_cluster_plot)
        hold on
        plot(proj_u(logical(in),1),proj_u(logical(in),2),'.r')
        hold off
            set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
            
       save('spikedata_removal.mat','in','-append');
end
    

% --- Executes during object creation, after setting all properties.
function popupmenu_selection_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_selection_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_select_spike.
function pushbutton_select_spike_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_spike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load('spikedata_removal.mat','proj_u','idx','axes_size','spikemat','in');

% if we==1
%     dm=[1 2];
% elseif we==2
%     dm=[1 3];
% elseif we==3
%     dm=[2 3];
% end
dm=[1 2];

if isempty(in)
    disp('Select Spikes Using Selection Method First')
    
else
 
    popup_label=[];
    cell_selectcluster = cellstr(get(handles.popupmenu_cluster_select,'string'));
    length_cell_selectcluster = length(cell_selectcluster);
    eval(sprintf('popup_label=''cluster%i'';',max(idx)+1))
    cell_selectcluster{length_cell_selectcluster + 1} = popup_label;
    set(handles.popupmenu_cluster_select, 'String', cell_selectcluster);
  
    idx(logical(in))=max(idx)+1;
    in_empty=idx==-1; % specifying the spike already removed by user
    
    
    axes(handles.axes_waveform_plot)
        
        if any(in_empty)
           plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
        else
           plot(spikemat,'color',[0.7 0.7 .7]) 
        end
            set(gca, 'XLim', [0 size(spikemat,1)]);
        
    axes(handles.axes_cluster_plot)
        
        if any(in_empty)
            plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
        else
            plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
        end
            set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
        
    in=[];
end

save('spikedata_removal.mat','idx','in','-append');


% --- Executes on selection change in popupmenu_cluster_select.
function popupmenu_cluster_select_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_cluster_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_cluster_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_cluster_select

load('spikedata_removal.mat','proj_u','idx','axes_size','spikemat');

% if we==1
%     dm=[1 2];
% elseif we==2
%     dm=[1 3];
% elseif we==3
%     dm=[2 3];
% end
dm=[1 2];

cell_selectcluster = cellstr(get(handles.popupmenu_cluster_select,'string'));
selected_cluster=get(handles.popupmenu_cluster_select, 'value');

set(handles.popupmenu_removal_method, 'value', 1);

in_empty=[]; in=[];

in_empty=idx==-1; % specifying the spike already removed by user

if selected_cluster==1
    
    axes(handles.axes_waveform_plot)
        
        if any(in_empty)
           plot(spikemat(:,~in_empty),'color',[0.7 0.7 0.7])
        else
           plot(spikemat,'color',[0.7 0.7 0.7]) 
        end
        set(gca, 'XLim', [0 size(spikemat,1)]);

    axes(handles.axes_cluster_plot)
        
        if any(in_empty)
            plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
        else
            plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
        end
            set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
            
else
    
popup_label=char(cell_selectcluster(selected_cluster));
cluster_no=str2num(popup_label( regexp(popup_label,'\d')));

in= idx==cluster_no;

  axes(handles.axes_waveform_plot)
  plot(spikemat(:,logical(in)),'r')
  set(gca, 'XLim', [0 size(spikemat,1)]);
  
  axes(handles.axes_cluster_plot)
  plot(proj_u(logical(in),1),proj_u(logical(in),2),'.r')
  set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]); 
  
end







% --- Executes during object creation, after setting all properties.
function popupmenu_cluster_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_cluster_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenu_removal_method.
function popupmenu_removal_method_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_removal_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_removal_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_removal_method

load('spikedata_removal.mat','proj_u','idx','axes_size','spikemat');

% if we==1
%     dm=[1 2];
% elseif we==2
%     dm=[1 3];
% elseif we==3
%     dm=[2 3];
% end
dm=[1 2];
selected_method=get(handles.popupmenu_removal_method, 'value');
cell_selectcluster = cellstr(get(handles.popupmenu_cluster_select,'string'));
selected_cluster=get(handles.popupmenu_cluster_select, 'value');

in_empty=idx==-1; % specifying the spike already removed by user

if selected_method==1 % No Selection Method Selected
    
    if selected_cluster==1 % No Cluster Selected
        
        axes(handles.axes_waveform_plot)
        
        if any(in_empty)
            plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
        else
            plot(spikemat,'color',[0.7 0.7 .7])
        end
        set(gca, 'XLim', [0 size(spikemat,1)]);
        
        axes(handles.axes_cluster_plot)
        
        if any(in_empty)
            plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
        else
            plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
        end
            set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
        
    else    % Cluster Selected
        
        in=[];
        popup_label=char(cell_selectcluster(selected_cluster));
        cluster_no=str2num(popup_label( regexp(popup_label,'\d')));
        
        in= idx==cluster_no;
        
        if any(in)==0
            
            cell_selectcluster(selected_cluster) =[];
            set(handles.popupmenu_cluster_select, 'String', cell_selectcluster);
            set(handles.popupmenu_cluster_select, 'value',1);
            set(handles.popupmenu_removal_method, 'value',1);
            
            axes(handles.axes_waveform_plot)
            
            if any(in_empty)
                plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
            else
                plot(spikemat,'color',[0.7 0.7 .7])
            end
            set(gca, 'XLim', [0 size(spikemat,1)]);
            
            axes(handles.axes_cluster_plot)
            
            if any(in_empty)
                plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
            else
                plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
            end
               set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
            
        else
            
            axes(handles.axes_waveform_plot)
            plot(spikemat(:,logical(in)),'k')
            set(gca, 'XLim', [0 size(spikemat,1)]);
            
            axes(handles.axes_cluster_plot)
            plot(proj_u(logical(in),1),proj_u(logical(in),2),'.k')
            set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
        end
        
    end
    
else   % Selection Method is Selected
    
    if selected_cluster ~= 1
        
        in=[]; in_cluster=[]; spikemat_cluster=zeros(size(spikemat)); proj_u_cluster = zeros(size(proj_u));
        popup_label=char(cell_selectcluster(selected_cluster));
        cluster_no=str2num(popup_label( regexp(popup_label,'\d')));
        in= idx==cluster_no;
        
        if any(in)==0   % No spike Present in the cluster
            
            cell_selectcluster(selected_cluster) =[];
            set(handles.popupmenu_cluster_select, 'String', cell_selectcluster);
            set(handles.popupmenu_cluster_select, 'value',1);
            set(handles.popupmenu_removal_method, 'value',1);
            
            axes(handles.axes_waveform_plot)
            
            if any(in_empty)
                plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
            else
                plot(spikemat,'color',[0.7 0.7 .7])
            end
            set(gca, 'XLim', [0 size(spikemat,1)]);
            
            axes(handles.axes_cluster_plot)
            
            if any(in_empty)
                plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
            else
                plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
            end
                set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
                
        else % spikes present in the cluster
            
            spikemat_cluster(:,in)=spikemat(:,in);
            proj_u_cluster(in,:)=proj_u(in,:);
            
            switch selected_method
                
                case 2
                    
                    axes(handles.axes_cluster_plot)
                    plot(proj_u(logical(in),1),proj_u(logical(in),2),'.k')
                    set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
                    
                    axes(handles.axes_waveform_plot)
                    plot(spikemat(:,logical(in)),'k')
                    set(gca, 'XLim', [0 size(spikemat,1)]);

                    X=[]; Y=[]; selected_points=[];  spikes_idx=[];
                    axes(handles.axes_cluster_plot)
                    [X,Y]= getpts;
                    selected_points=[X,Y];
                    spikes_idx = knnsearch(proj_u_cluster,selected_points);
                    
                    in_cluster=zeros(1,size(spikemat,2));
                    in_cluster(spikes_idx)=1;
                    
                    axes(handles.axes_cluster_plot)
                    hold on
                    plot(proj_u_cluster(spikes_idx,1),proj_u_cluster(spikes_idx,2),'*')
                    hold off
                    set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
                    
                    axes(handles.axes_waveform_plot)
                    hold on
                    plot(spikemat_cluster(:,spikes_idx))
                    hold off
                    set(gca, 'XLim', [0 size(spikemat,1)]);

                    save('spikedata_removal.mat','in_cluster','-append');
                    
                    
                case 3
                    
                    axes(handles.axes_waveform_plot)
                    plot(spikemat(:,logical(in)),'k')
                    set(gca, 'XLim', [0 size(spikemat,1)]);

                    axes(handles.axes_cluster_plot)
                    plot(proj_u(logical(in),1),proj_u(logical(in),2),'.k')
                    set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
                    
                    X=[]; Y=[]; spikes_idx=[];
                    axes(handles.axes_waveform_plot)
                    [X,Y]= getpts;
                    
                    for k=1:length(X)
                        part_data=[];
                        
                        if floor(X(k))~=X(k) && floor(X(k))~=0
                            part_data=spikemat_cluster(floor(X(k)):ceil(X(k)),:);
                        elseif floor(X(k))== 0
                            part_data=spikemat_cluster(1,:);
                        else
                            part_data=spikemat_cluster(X(k),:);
                        end
                        
                        
                        spikes_idx(k) = knnsearch(part_data',repmat(Y(k),1,size(part_data,1)));
                        
                    end
                    
                    in_cluster=zeros(1,size(spikemat,2));
                    in_cluster(spikes_idx)=1;
                    
                    axes(handles.axes_waveform_plot)
                    hold on
                    plot(spikemat_cluster(:,spikes_idx))
                    hold off
                    set(gca, 'XLim', [0 size(spikemat,1)]);
                    
                    axes(handles.axes_cluster_plot)
                    hold on
                    plot(proj_u_cluster(spikes_idx,1),proj_u_cluster(spikes_idx,2),'*')
                    hold off
                    set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
                    
                    save('spikedata_removal.mat','in_cluster','-append');
                    
                    
                case 4
                    
                    rect=[]; part_data=[];  in_cluster=[];
                    
                    axes(handles.axes_waveform_plot)
                    plot(spikemat(:,logical(in)),'k')
                    set(gca, 'XLim', [0 size(spikemat,1)]);
                    
                    axes(handles.axes_cluster_plot)
                    plot(proj_u(logical(in),1),proj_u(logical(in),2),'.k')
                    set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
                    
                    rect=getrect(handles.axes_waveform_plot);
                    if rect(1)>=1
                        part_data=spikemat_cluster(floor(rect(1)):ceil(rect(1)+rect(3)),:);
                    else
                        part_data=spikemat_cluster(1:ceil(rect(1)+rect(3)),:);
                    end
                    for i=1:size(part_data,2)
                        in_cluster(i)=any(part_data(:,i)>=rect(2)) & any(part_data(:,i)<=rect(2)+rect(4));
                    end
                    
                    axes(handles.axes_waveform_plot)
                    hold on
                    plot(spikemat_cluster(:,logical(in_cluster)),'r')
                    hold off
                    set(gca, 'XLim', [0 size(spikemat,1)]);
                    
                    axes(handles.axes_cluster_plot)
                    hold on
                    plot(proj_u_cluster(logical(in_cluster),1),proj_u_cluster(logical(in_cluster),2),'*r')
                    hold off
                    set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
                    
                    save('spikedata_removal.mat','in_cluster','-append');
                    
                    
                case 5
                    
                    in_cluster=[];
                    axes(handles.axes_cluster_plot)
                    plot(proj_u(logical(in),1),proj_u(logical(in),2),'.k')
                    set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
                    
                    axes(handles.axes_waveform_plot)
                    plot(spikemat(:,logical(in)),'k')
                    set(gca, 'XLim', [0 size(spikemat,1)]);
                    
                    axes(handles.axes_cluster_plot)
                    convex_boundary_x=[]; convex_boundary_y=[]; k=[];
                    [convex_boundary_x,convex_boundary_y] = ginput();
                    convex_boundary_x = [convex_boundary_x;convex_boundary_x(1)];
                    convex_boundary_y = [convex_boundary_y;convex_boundary_y(1)];
                    k = convhull(convex_boundary_x,convex_boundary_y);
                    convex_boundary_x = convex_boundary_x(k);
                    convex_boundary_y = convex_boundary_y(k);
                    in_cluster = inpolygon(proj_u_cluster(:,1),proj_u_cluster(:,2),convex_boundary_x,convex_boundary_y);
                    
                    axes(handles.axes_waveform_plot)
                    hold on
                    plot(spikemat_cluster(:,logical(in_cluster)),'r')
                    hold off
                    set(gca, 'XLim', [0 size(spikemat,1)]);
                    
                    axes(handles.axes_cluster_plot)
                    hold on
                    plot(proj_u_cluster(logical(in_cluster),1),proj_u_cluster(logical(in_cluster),2),'*r')
                    hold off
                    set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
                    
                    save('spikedata_removal.mat','in_cluster','-append');
            end
            
        end
        
    else % Selection Method Selected Without Selecting Cluster
        
        disp('Select a Cluster From Cluster List Before Selecting Removal Method')
        
    end
end




% --- Executes during object creation, after setting all properties.
function popupmenu_removal_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_removal_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in pushbutton_remove_from_cluster.

function pushbutton_remove_from_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove_from_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushbutton_remove_from_cluster


load('spikedata_removal.mat','proj_u','idx','axes_size','spikemat','in_cluster');

% if we==1
%     dm=[1 2];
% elseif we==2
%     dm=[1 3];
% elseif we==3
%     dm=[2 3];
% end
dm=[1 2];
selected_method=get(handles.popupmenu_removal_method, 'value');
cell_selectcluster = cellstr(get(handles.popupmenu_cluster_select,'string'));
selected_cluster=get(handles.popupmenu_cluster_select, 'value');

in_empty=idx==-1; % specifying the spike already removed by user


if selected_cluster~=1 && selected_method~=1
    
    if any(in_cluster)~=1
        
        disp('No Spike Selected From the Cluster for Removal.....Select a Spike from the cluster with a Selection Method')
        
    else
        
        popup_label=char(cell_selectcluster(selected_cluster));
        cluster_no=str2num(popup_label( regexp(popup_label,'\d')));
        
        idx(find(in_cluster==1))=0;
        save('spikedata_removal.mat','idx','-append');
        in= idx==cluster_no;
        
        if any(in)==0   % No spike Present in the cluster
            
            cell_selectcluster(selected_cluster) =[];
            set(handles.popupmenu_cluster_select, 'String', cell_selectcluster);
            set(handles.popupmenu_cluster_select, 'value',1);
            set(handles.popupmenu_removal_method, 'value',1);
            
            axes(handles.axes_waveform_plot)
            
            if any(in_empty)
                plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
            else
                plot(spikemat,'color',[0.7 0.7 .7])
            end
            set(gca, 'XLim', [0 size(spikemat,1)]);
            
            axes(handles.axes_cluster_plot)
            
            if any(in_empty)
                plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
            else
                plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
            end
            set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
            
        else
            
            axes(handles.axes_waveform_plot)
            plot(spikemat(:,logical(in)),'color',[0.7 0.7 .7])
            set(gca, 'XLim', [0 size(spikemat,1)]);
            
            axes(handles.axes_cluster_plot)
            plot(proj_u(logical(in),1),proj_u(logical(in),2),'.k')
            set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
            
            
            set(handles.popupmenu_removal_method, 'value', 1);
            
        end
        
    end
    
else
    
    disp('Select Either Cluster Or Method of Removal Or Both Before Selecting This....')
    
end





% --- Executes on button press in pushbutton_remove.

function pushbutton_remove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load('spikedata_removal.mat','proj_u','idx','spikemat','axes_size');

% if we==1
%     dm=[1 2];
% elseif we==2
%     dm=[1 3];
% elseif we==3
%     dm=[2 3];
% end
dm=[1 2];
cell_selectcluster = cellstr(get(handles.popupmenu_cluster_select,'string'));
selected_cluster=get(handles.popupmenu_cluster_select, 'value');

if selected_cluster~=1 
    
    popup_label=char(cell_selectcluster(selected_cluster));
    cluster_no=str2num(popup_label( regexp(popup_label,'\d')));

    in= idx==cluster_no;
    idx(find(in==1))=-1;
    spikemat(:,in)=spikemat(:,in)==0;
    proj_u(in,:)=proj_u(in,:)==0;
    
    save('spikedata_removal.mat','idx','spikemat','proj_u','-append'); 
     
    cell_selectcluster(selected_cluster) =[];
    set(handles.popupmenu_cluster_select, 'String', cell_selectcluster);
    set(handles.popupmenu_cluster_select, 'value',1);
    set(handles.popupmenu_removal_method, 'value',1);
    set(handles.popupmenu_selection_method, 'value',1);
            
    in_empty=idx==-1; % specifying the spike already removed by user
    
    axes(handles.axes_waveform_plot)
    plot(spikemat(:,~in_empty),'color',[0.7 0.7 .7])
    set(gca, 'XLim', [0 size(spikemat,1)]);

    axes(handles.axes_cluster_plot)
    plot(proj_u(~in_empty,1),proj_u(~in_empty,2),'.k')
    set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);
  
else
    
    disp('No Cluster Selected : Select Cluster From Cluster List Before Selecting This Button......')
    
end    
    
     
    

% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load('spikedata_removal.mat','spikemat','spkid','spktime','idx');

% if we==1
%     dm=[1 2];
% elseif we==2
%     dm=[1 3];
% elseif we==3
%     dm=[2 3];
% end
dm=[1 2];
idx_index=find(idx==-1);

idx=zeros(size(idx));

if ~isempty(idx_index)
    idx(idx_index)=-1; 
end

% save('spikedata.mat','idx','-append'); 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(fields(spkid))
    
spk_field=[]; spk_time=[];

eval(sprintf('spk_field=spkid.iter%i;',i))
eval(sprintf('spk_time=spktime.iter%i;',i))
eval(sprintf('spkid.iter%i=[];',i))
eval(sprintf('spktime.iter%i=[];',i))

if any(ismember(spk_field,idx_index))
    
    spk_time(ismember(spk_field,idx_index))=[];
    spk_field(ismember(spk_field,idx_index))=[];
    
end

eval(sprintf('spkid.iter%i=spk_field;',i))
eval(sprintf('spktime.iter%i=spk_time;',i))

end

start=1;

for i=1:length(fields(spkid))
    
  field_len=0;  

eval(sprintf('field_len=length(spkid.iter%i);',i))
eval(sprintf('spkid.iter%i=[];',i))
eval(sprintf('spkid.iter%i=start:1:start + field_len-1;',i))

start=start+field_len;

end

spikemat(:,idx_index)=[];


handles.spikemat=spikemat;
guidata(hObject, handles);
handles.spkid=spkid;
guidata(hObject, handles);
handles.spktime=spktime;
guidata(hObject, handles);
handles.idx=idx; 
guidata(hObject, handles);


close(handles.spike_removal_gui);


% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


load('spikedata_removal.mat','spikemat_copy','proj_u_copy','axes_size');

% if we==1
%     dm=[1 2];
% elseif we==2
%     dm=[1 3];
% elseif we==3
%     dm=[2 3];
% end
dm=[1 2];
spikemat=spikemat_copy;
proj_u=proj_u_copy;
idx=zeros(size(spikemat,2),1);
clusters=[];


    set(handles.popupmenu_cluster_select, 'String', 'Select Cluster To Remove');
    set(handles.popupmenu_cluster_select, 'value',1);
    set(handles.popupmenu_removal_method, 'value',1);
    set(handles.popupmenu_selection_method, 'value',1);

    save('spikedata_removal.mat','proj_u','idx','spikemat','-append');

axes(handles.axes_cluster_plot)
plot(proj_u(:,dm(1)),proj_u(:,dm(2)),'.k')
set(gca, 'XLim', [axes_size(1) axes_size(3)], 'YLim', [axes_size(2) axes_size(4)]);

axes(handle.axes_waveform_plot)
plot(spikemat,'color',[0.7 0.7 0.7])
set(gca, 'XLim', [0 size(spikemat,1)]);


% --- Executes when user attempts to close spike_removal_gui.
function spike_removal_gui_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to spike_removal_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    
% Hint: delete(hObject) closes the figure
delete(hObject);

end






