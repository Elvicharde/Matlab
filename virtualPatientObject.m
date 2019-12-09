classdef virtualPatientObject < handle
    % This is a virtual patient object.
    %
    % The virtual patient design is based on the Hovorka model. An instance
    % of this class has options for including a controller and exercise for
    % simulation i.e. incorporation of exercise dynamics.
    %
    % To create a virtual patient, use the following syntax:
    %
    %  objectName = virtualPatientObject('name', age, 'sex', BW)
    %      where: The simulated patient's name = 'name'
    %             age is a numeric value
    %             sex is either 'Male' or 'Female'
    %             BW represents the body weight in Kg.
    % 
    %   Example: Patient1 = virtualPatientObject('Jin', 50, 'Male', 77)
    %
    % This assigns the object 'Patient 1' with attributes specified within
    % the virtualPatientObject's arguments. 
    %
 %%   
    properties 
    % Here is a declaration of all the Public properties of the
    % virtual patient.
        
        % Bio data properties
        name
        age 
        sex
        body_Weight 
        diabeties_Variant = 'Type 1 Diabetes Mellitus'
        
        % Other Properties      
        insulin_Tolerance     
        controller_Status = 'Detached'  % Default state
        exercise = 'False'              % Default state
    end
   
    properties 
        foodIn                 % Meal input quantity [grammes]
        insulin_In             % Insulin rate from controller [mU/min].
    end
    
    properties (Access = protected)
    % These are the properties that are internal to this object. They are 
    % used for calculating the values of some of the public properties.
    
    % Food compartment properties       
        Dg                     % Glucose equivalent of ingested 
                               ... carbohydrate [mmol/min].
        D1                     % First glucose absorption compartment
        D2                     % second glucose absorption compartment
        Ug                     % output to bloodstream   
        
    % Insulin compartment properties          
        S1                     % First Insulin absorption compartment [mU]
        S2                     % Second Insulin absorption compartment [mU]
        UI                     % Absorption rate to bloodstream [mU/min]
    end
    
    properties (Access = protected)
    % Physical attributes/states  
        Q1                     % Glucose in accessible plasma [mmol]
        Q2                     % Glucose in non-accessible 
                               ...compartment [mmol]   
    end
    
    properties (Constant)
    % These are already defined values that do not change.    
        Mwg  = 180.16;         % Molecular weight of glucose.
        Td   = 40;             % Time-to-maximum of carbohydrate
                               ... absorption (mins)
        Ag   = 0.8;            % Glucose utilization factor
        Ts   = 55;             % Time-to-maximum of insulin 
                                ... absorption & action (mins)
        Ke   =  0.138          % Elimination rate of insulin [min-1]
        k12  =  0.066          % Fractional transfer rate Q1-Q2 [min-1]
        ka1  =  0.006          % Deactivation rate parameter [min-1]
        ka2  =  0.06           % Deactivation rate parameter [min-1]
        ka3  =  0.03           % Deactivation rate parameter [min-1]
        Si1  =  0.00512        % Transport Insulin sensitivity [min-1/mU/L]
        Si2  =  0.00082        % Disposal Insulin sensitivity [min-1/mU/L]
        Si3  =  0.052          % EGP insulin activity [L/mU]
        Vg   =  0.16           % Distribution volume of glucose [L]
        Vi   =  0.12           % Distribution volume of insulin [L]
        Fo1  =  0.0111         % Glucose uptake by CNS  
                               ... (Central nervous system)[mmol/min]
        EGPo =  0.0161         % Endogenous glucose production (Hepatic)
                               ...[mmol/min] 
    end
    
     properties
        % Parameters 
        Fo1c                   % Variable Fo1
        Fr                     % Renal glucose uptake dispense [mmol/L]
        
    end
    
    properties (Dependent)
    % These properties vary based on the value of other properties. 
    % Some are also event triggers.
      
    % Parameters
        plasma_InsulinConcentration     % I(t) [mU/L]
        plasma_GlucoseConcentration     % G(t) [mmol/L]        
        EGP                             % Hepatic endogenous glucose   
                                            
    % Insulin-on-glucose effect parameters                                    
        x1                              % Glucose distribution factor
        x2                              % Glucose disposal factor
        x3                              % Glucose production factor
    end
    
    %%    
    events
    % These are occurences that require an action: mostly control actions
        high_Blood_Sugar
        low_Blood_Sugar
        high_Insulin_Level
        controller_Action
    end
    
    %% Here goes the static methods for all instances of this object.
    methods (Static)
        % Adding the object constructor
        function obj = virtualPatientObject(varargin)
            % Class constructor: This is called at any instance of the
            ... class.
                
        if isempty(varargin)
            disp((sprintf(['Virtual Patient created without details!\n\n',...
                'There are four property values: name,age,sex(Male/',...
                'Female),and body_Weight.\nUse the syntax objname.property',...
                ' method to assign values.\n'])));            
        elseif nargin
            if nargin < 4
                disp([newline 'One or more parameters missing' newline])
                clear 
                
            elseif nargin > 4
                disp([newline 'Too many inputs' newline])
                clear 
            else
                obj.name = lower(varargin{1});
                obj.age = double(varargin{2});
                obj.sex = varargin{3};
                obj.body_Weight = double(varargin{4});
                
            end
        end
        end
    end
    
%    debug status: working!
%% Here goes the dynamic methods for all instances of this object.
    
    methods     
%%      Some 'Get'-er methods  
%--------------------------------------------------------------------------       
        function out = get.plasma_InsulinConcentration(obj)
        % This method finds the plasma insulin concentration of this
        % virtual patient during the simulation. It updates the value of 
        % the plasma_InsulinConcentration property, which is then displayed
        
        % Setting condition for plasma insulin call.
        infusion_time = 10; % 10 minutes for constant drug infusion
        Ui = [0,0];
        
        if ~isempty(obj.insulin_In)
            for i = 1:300
                if i <= infusion_time
                    % Ui  [Mu/min]
                    Ui(i) = (obj.insulin_In*1000)/infusion_time;
                else
                    Ui(i) = 0;
                end
            end
            
            Ui = [0,Ui];
            Ut = linspace(0,300,301);           % Interpolation of U.            
            
        else
            for i = 1:300
                Ui(i) = 0;
            end
            Ui = [0,Ui];
            Ut = linspace(0,300,301);           % Interpolation of D.
        end
        absorbInsulin(obj,Ui,Ut);          % Calling the absorb method
        
%                        Computation of plasma insulin        

        len = 301 - length(obj.UI);
        pad = zeros(len,1);
        flowrate = [obj.UI;pad]';
        Ut = linspace(0,300,301);       % Vector for interpolation.

        % Setting the time range:
        t_start = 0;
        t_final = 300;

        % Initial condition(s)
        I_init = 5.8;   % basal insulin value [mU/L].

        % Solving the ODE using Dormand-Prince 45 numerical solver

        [T, I] = ode45(@(t,i) plasma_conc(obj,t,i,Ut,flowrate)...
            ,[t_start, t_final],I_init);

        
        out = I;
        
        % Plots
        figure()
        plot(T,I,'g', 'LineWidth',0.8)
%         ylabel('I(t) [mU/L]','FontSize',12)
%         xlabel('Time (hrs)','FontAngle','italic')
%         xticks(0:30:400)
%         xticklabels({'8:00 AM','8:30 AM','9:00 AM','9:30 AM'...
%             '10:00 AM','10:30 AM','11:00 AM','11:30 AM','12:00 PM', ...
%             '12:30 PM', '1:00 PM','1:30 PM'})
%         title ('Plasma Insulin concentration')
%         grid on

%        Differential equation function for I(t)

                    function didt = plasma_conc(obj,t,i,Ut,flowrate)
                        VI = obj.Vi*obj.body_Weight; 
                        U_int = interp1(Ut,flowrate,t); % interpolating 
                        didt = (U_int/VI) - (obj.Ke*i);
                    end

        end
% debug status: working!        
%--------------------------------------------------------------------------       
        function out = get.plasma_GlucoseConcentration(obj)
        % This method finds the plasma glucose concentration of this
        % virtual patient during the simulation. It updates the value of 
        % the plasma_GlucoseConcentration property, which is then displayed
        
        % declaring variables and dependencies
        VG = obj.Vg * obj.body_Weight;
                
        % This method gets the amount of glucose in the plasma.  
        
        % Preparing all vectors for Interpolation 
        
        X1 = obj.x1;
        len1 = length(X1);
        Ix1 = linspace(0,len1-1,len1);
        
        X2 = obj.x2;
        len2 = length(X2);
        Ix2 = linspace(0,len2-1,len2);
        
        if isempty(obj.foodIn) || obj.foodIn == 0
            ug = zeros(1,300)';
        else
            ug = obj.Ug;
        end
        
        len3 = length(ug);
        len = 301 - len3;
        pad = zeros(len,1);
        ug  = [ug;pad]';
        IUg = linspace(0,300,301);              % Vector for interpolation.
        
        eGP = obj.EGP;
        len4 = length(eGP);
        IEGP = linspace(0,len4-1,len4);
        
        
        
        % Now the solution
        
        t_start = 0;
        t_final = 300;
        
        % Initial condition(s):
        Q1_init = 54;   %[mmol]
        Q2_init = 0;    %[mmol]
        
        [T, Q] = ode45(@(t,q) gluc(obj,t,q,Ix1,X1,Ix2,X2,IUg,ug,IEGP,...
                 eGP),[t_start, t_final],[Q1_init,Q2_init]);
        
        obj.Q1 = Q(:,1);
        obj.Q2 = Q(:,2);
        G2     = obj.Q1/VG;
        out    = G2;
        
        % Call plot
        
        x = T;
        figure()
        plot(x,G2,'g', 'LineWidth',0.8)
        ylabel('G [mmol/L]','FontWeight','Bold','FontSize',14)
        xlabel('Time (hrs)','FontAngle','Italic')
        xticks(0:30:400)
        xticklabels({'8:00 AM','8:30 AM','9:00 AM','9:30 AM'...
            '10:00 AM','10:30 AM','11:00 AM','11:30 AM','12:00 PM', ...
            '12:30 PM', '1:00 PM','1:30 PM'})
        title ('Plasma glucose concentration')
        grid on
        
        figure()
        clf;
        subplot(2, 1, 1)
        plot (x,obj.Q1,'r','LineWidth',0.8)
        ylabel(' Q_1 [mmol]','FontSize',12)
        xlabel('Time (hrs)')
        xticks(0:30:400)
        xticklabels({'8:00 AM','8:30 AM','9:00 AM','9:30 AM'...
            '10:00 AM','10:30 AM','11:00 AM','11:30 AM','12:00 PM','12:30 PM',...
            '1:00 PM','1:30 PM'})
        title (['Plasma Glucose for ',num2str(obj.foodIn), 'g of CHO'])
        grid on
        hold on
        
        subplot(2, 1, 2)
        plot (x,obj.Q2,'LineWidth',0.8)
        ylabel('Q_2 [mmol]','FontSize',12)
        xlabel('Time (hrs)')
        xticks(0:30:400)
        xticklabels({'8:00 AM','8:30 AM','9:00 AM','9:30 AM'...
            '10:00 AM','10:30 AM','11:00 AM','11:30 AM','12:00 PM','12:30 PM',...
            '1:00 PM','1:30 PM'})
        grid on
        hold on
      
        % Actual function for solving the ODE
         
            function dqdt = gluc(obj,t,q,Ix1,X1,Ix2,X2,IUg,ug,IEGP,eGP)
                
                x1_int  = interp1(Ix1,X1,t);   % interpolating the  
                                              ...input data for evaluation.
                x2_int  = interp1(Ix2,X2,t); 
                Ug_int  = interp1(IUg,ug,t); 
                EGP_int = interp1(IEGP,eGP,t); 
                
                % The conditional parameters
                G = q(1)/(obj.Vg *obj.body_Weight);
                
                
                % Fo1c
                if G >= 4.5 % 4.5 mmol/L
                    fo1c     = obj.Fo1;                    
                else
                    fo1c = obj.Fo1*(G/4.5);                    
                end
                
                % FR
                if G >= 9 % 4.5 mmol/L
                    fR = 3e-3*(G - 9)*obj.Vg;                    
                else
                    fR = 0;                    
                end
                
                % Assign some more variables for each compartment
                
                q1 = q(1);
                q2 = q(2);
                
                
                % for Q1
                % The excreted elements/out-flow
                out_flow = fo1c + fR + (x1_int * q2);
                
                % The in-flow
                in_flow = Ug_int + EGP_int + (obj.k12 * q2);
                
                % for Q2
                % The excreted elements/out-flow
                out_flow2 = (obj.k12 + x2_int)*q2;
                
                % The in-flow
                in_flow2 = x1_int * q1;
                
                
        % The output from the ODE function is a COLUMN vector with 2 rows
                
                dqdt    = zeros(2,1);
                dqdt(1) = in_flow - out_flow;  
                dqdt(2) = in_flow2 - out_flow2;
            end      
        
        end
%--------------------------------------------------------------------------        
        function out = get.EGP(obj)
        % This gets the current value of endogenous glucose production 
        % activity of the liver. It depends on the EGPo parameter and 
        % the endogenous glucose production parameter 'x3'.    
        
        val = [0,0];
        for i = 1:length(obj.x3)
            val(i) = obj.EGPo* (1 - obj.x3(i));
        end
        out = val';
 
        end 
%--------------------------------------------------------------------------        
        function out = get.x1(obj)
        % This gets the current value of the glucose distribution parameter
        % It depends on the plasma insulin concentration and has an effect
        % on the overall plasma glucose concentration.
        I = obj.plasma_InsulinConcentration;
        len = length(I);
        It = linspace(0,len-1,len);
        
        % Setting the time range:
        t_start = 0;
        t_final = len-1;
        
        % Initial condition(s):
        x1_init = 0.0295;
        
        [T, X1] = ode45(@(t,x1) x_val(obj,t,x1,It,I)...
            ,[t_start, t_final],x1_init);
        
        out = X1;
        
        %  differential equation function for I(t)
            function dxdt = x_val(obj,t,x1,It,I)
                I_int = interp1(It,I,t); % interpolating 
                kb1 = obj.ka1 * obj.Si1;
                dxdt = -obj.ka1*x1 + kb1*I_int;
            end
        end 
        
%--------------------------------------------------------------------------        
        function out = get.x2(obj)
        % This gets the current value of the glucose disposal parameter
        % It depends on the plasma insulin concentration and has an effect
        % on the overall plasma glucose concentration.
        
        I = obj.plasma_InsulinConcentration;
        len = length(I);
        It = linspace(0,len-1,len);
        
        % Setting the time range:
        t_start = 0;
        t_final = len-1;
        
        % Initial condition(s):
        x2_init = 0.0047; 
        
        [T, X2] = ode45(@(t,x2) x_val(obj,t,x2,It,I)...
            ,[t_start, t_final],x2_init);
        
        out = X2;
        
        %  differential equation function for I(t)
            function dxdt = x_val(obj,t,x2,It,I)
                I_int = interp1(It,I,t); % interpolating 
                kb2 = obj.ka2 * obj.Si2;
                dxdt = -obj.ka2*x2 + kb2*I_int;
            end
          
        end 
%--------------------------------------------------------------------------        
        function out = get.x3(obj)
        % This gets the current value of the endogenous glucose production
        % parameter. It depends on the plasma insulin concentration and has
        % an effect on the overall plasma glucose concentration.    
            
        I = obj.plasma_InsulinConcentration;
        len = length(I);
        It = linspace(0,len-1,len);
        
        % Setting the time range:
        t_start = 0;
        t_final = len-1;
        
        % Initial condition(s):
        x3_init = 0.2996; 
        
        [T, X3] = ode45(@(t,x3) x_val(obj,t,x3,It,I)...
            ,[t_start, t_final],x3_init);
        
        out = X3;
        
        %  differential equation function for I(t)
            function dxdt = x_val(obj,t,x3,It,I)
                I_int = interp1(It,I,t); % interpolating 
                kb3 = obj.ka1 * obj.Si3;
                dxdt = -obj.ka3*x3 + kb3*I_int;
            end
            
        end 
     
%% Some Activity methods
%--------------------------------------------------------------------------
% This function runs a simulation of the object 
        function out = simulate(obj,controller_Status,exercise, sim_Time)
        % This method simulates the virtual patient object with
        % options for an attached controller or exercise activity. 
            
            
            
        end        

%--------------------------------------------------------------------------
        function eat(obj,meal)
        % This method calls the object to eat a meal i.e. glucose
        % disturbance. It calculates the glucose value of the meal and then
        % passes it to the glucose absorption compartment.
        
        if strcmpi(meal,'Breakfast')
            obj.foodIn = 25;   % 25.0 grammes of carbohydrate from
                               % eating a Sandwich and Raisins.
            duration  = 10;    % Time to complete meal        
        elseif strcmpi(meal,'Lunch')
            obj.foodIn = 57;   % 57.0 grammes of carbohydrate from eating 
                               % Rice, cooked carrots, honey dew, and curry 
                               % korma.
            duration  = 20;    % Time to complete meal            
        elseif strcmpi(meal,'Dinner')
            obj.foodIn = 86.5; % 86.5 grammes of carbohydrate eating Rice,
                               % fish fingers, cooked cabbage, and banana.
            duration  = 30;    % Time to complete meal                    
        end
        
        for i = 1:300
            if i <= duration
                % Dg  [g/min]
                obj.Dg(i) = (obj.foodIn*1000)/(duration*obj.Mwg);
            else
                obj.Dg(i) = 0;
            end
        end 
        
        Di = [0,obj.Dg];
        Dt = linspace(0,300,301);     % Vectors for interpolation of D.
        
        absorbMeal(obj,Di,Dt);
                       
        end
% Debug status: This works!
%--------------------------------------------------------------------------
%      
%         function out = exercise()
%         % This method calls the object to eat a meal i.e. glucose
%         % disturbance. It calculates the glucose value of the meal and then
%         % passes it to the glucose absorption compartment.
%         
%         
%         
%         end
%%      Some Absorption methods
%--------------------------------------------------------------------------
        function absorbMeal(obj,Di,Dt)
        % This is the object meal absorption method. It also has an 
        % overload of the plot function to display meal profile.
        % 
        % Parameter call is from within either the eat() method.
        
        % Setting the time range:
        t_start = 0;
        t_final = 300;
        
        % Initial condition(s):
        D1_init = 0;
        D2_init = 0;
        
        [T, D] = ode45(@(t,d) food(obj,t,d,Dt,Di),[t_start, t_final]...
            ,[D1_init, D2_init]);
        
        obj.D1 = D(:,1);
        obj.D2 = D(:,2);
        obj.Ug = obj.D2/obj.Td;
        
        % Call plot
        
        x = T;
        figure()
        clf;
        subplot(2, 1, 1)
        plot (x,obj.D1,'r','LineWidth',0.8)
        ylabel(' D1 [mmol]','FontSize',12)
        xlabel('Time (hrs)')
        xticks(0:30:400)
        xticklabels({'8:00 AM','8:30 AM','9:00 AM','9:30 AM'...
            '10:00 AM','10:30 AM','11:00 AM','11:30 AM','12:00 PM','12:30 PM',...
            '1:00 PM','1:30 PM'})
        title (['Glucose absorption profile for ',num2str(obj.foodIn), 'g of CHO'])
        grid on
        hold on
        
        subplot(2, 1, 2)
        plot (x,obj.D2,'LineWidth',0.8)
        ylabel('D2 [mmol]','FontSize',12)
        xlabel('Time (hrs)')
        xticks(0:30:400)
        xticklabels({'8:00 AM','8:30 AM','9:00 AM','9:30 AM'...
            '10:00 AM','10:30 AM','11:00 AM','11:30 AM','12:00 PM','12:30 PM',...
            '1:00 PM','1:30 PM'})
        grid on
        hold on
        
        figure()
        plot(x,obj.Ug,'g', 'LineWidth',0.8)
        ylabel('Ug [mmol/min]','FontWeight','Bold','FontSize',14)
        xlabel('Time (hrs)','FontAngle','Italic')
        xticks(0:30:400)
        xticklabels({'8:00 AM','8:30 AM','9:00 AM','9:30 AM'...
            '10:00 AM','10:30 AM','11:00 AM','11:30 AM','12:00 PM', ...
            '12:30 PM', '1:00 PM','1:30 PM'})
        title ('Glucose gut absorption rate')
        grid on
      
        % Actual function for solving the ODE
         
            function dddt = food(obj,t,d,Dt,Di)
                
                D_int = interp1(Dt,Di,t); % interpolating the input data 
                                          ...for evaluation.
                
        % Assign some more variables for each compartment
                
                d1 = d(1);
                d2 = d(2);
                
        % The output from the ODE function is a COLUMN vector with 2 rows
                
                dddt    = zeros(2,1);
                dddt(1) = (D_int*obj.Ag) - (d1/obj.Td);
                dddt(2) = (d1 - d2)/obj.Td;
                
            end      
        
        end
        
% Debug status: This works!
%-------------------------------------------------------------------------- 
        function  out = absorbInsulin(obj,Ui,Ut)
        % This function is called to absorb the insulin dosage prescribed 
        % by the attached controller. If no controller is attached, it does
        % not have a value. Since a simulation could be run without insulin
        % the plot function here is set to only display when there are
        % values to display.
        
        % Setting the time range:
        t_start = 0;
        t_final = 300;
        
        % Initial condition(s):
        U1_init = 0;
        U2_init = 0;
        
        % Solving the ODE(s) using Dormand-Prince 45 numerical solver
        
        [T, S] = ode45(@(t,s) insulin_action(obj,t,s,Ut,Ui),[t_start,... 
                                           t_final],[U1_init, U2_init]);
        
        obj.S1 = S(:,1);
        obj.S2 = S(:,2);
        obj.UI = obj.S2/obj.Ts;
        out = obj.UI;
        
%         % Plots
%         if ~isempty(obj.insulin_In)
%             x = T;
%             figure()
%             clf;
%             subplot(2, 1, 1)
%             plot (x,obj.S1,'r','LineWidth',0.8)
%             ylabel(' U1(t) [mU]','FontWeight','Bold','FontSize',14)
%             xlabel('Time (hrs)')
%             xticks(0:30:400)
%             xticklabels({'8:00 AM','8:30 AM','9:00 AM','9:30 AM'...
%                 '10:00 AM','10:30 AM','11:00 AM','11:30 AM','12:00 PM',...
%                 '12:30 PM','1:00 PM','1:30 PM'})
%             title (['Absorption profile for ',num2str(obj.insulin_In),...
%                                                       'U of Insulin'])
%             grid on
%             hold on
% 
%             subplot(2, 1, 2)
%             plot (x,obj.S2,'LineWidth',0.8)
%             ylabel('U2(t) [mU]','FontWeight','Bold','FontSize',14)
%             xlabel('Time (hrs)')
%             xticks(0:30:400)
%             xticklabels({'8:00 AM','8:30 AM','9:00 AM','9:30 AM'...
%                 '10:00 AM','10:30 AM','11:00 AM','11:30 AM','12:00 PM',...
%                 '12:30 PM','1:00 PM','1:30 PM'})
%             grid on
%             hold on
% 
%             figure()
%             plot(x,obj.UI,'g', 'LineWidth',0.8)
%             ylabel('Ug(t) [mU/min]','FontSize',12)
%             xlabel('Time (hrs)','FontAngle','italic')
%             xticks(0:30:400)
%             xticklabels({'8:00 AM','8:30 AM','9:00 AM','9:30 AM'...
%                 '10:00 AM','10:30 AM','11:00 AM','11:30 AM','12:00 PM',...
%                 '12:30 PM', '1:00 PM','1:30 PM'})
%             title ('Insulin gut absorption rate')
%             grid on
%         end
       
        % differential equation function for S1 & S2
            function dsdt = insulin_action(obj,t,s,Ut,Ui)
                
                U_int = interp1(Ut,Ui,t); % interpolating the input data 
                
         % Assign some more variables for each compartment
                
                s1 = s(1);
                s2 = s(2);
                
          % The output from the ODE function is a COLUMN vector with 2 rows
                
                dsdt    = zeros(2,1);
                dsdt(1) = U_int - (s1/obj.Ts);
                dsdt(2) = (s1 - s2)/obj.Ts;
                
            end
        end
        
    end
%    debug status: working!    
end