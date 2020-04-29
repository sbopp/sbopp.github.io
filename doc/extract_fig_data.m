%%
% MATLAB Data Extraction from .fig Files
% Written by Steven E. Bopp April 28, 2020

%% SiO2 Au SiO2 Wells

h_SiO2_Au3nm_SiO2=openfig('/Users/mainuser/Desktop/Au SHG Paper Figures/SiO2_Au3nm_SiO2_with_fermi_energy.fig');
h_SiO2_Au3nm_SiO2=findobj(gca,'Type','Line');
x_SiO2_Au3nm_SiO2=get(h_SiO2_Au3nm_SiO2,'Xdata');
y_SiO2_Au3nm_SiO2=get(h_SiO2_Au3nm_SiO2,'Ydata');

%y_Fermi_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{1,1};
%y_Elvl10_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{2,1};
y_Elvl9_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{3,1};
y_Elvl8_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{4,1};
y_Elvl7_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{5,1};
%y_Elvl6_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{6,1};
%y_Elvl5_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{7,1};
%y_Elvl4_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{8,1};
%y_Elvl3_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{9,1};
%y_Elvl2_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{10,1};
%y_Elvl1_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{11,1};
%y_Well_SiO2_Au3nm_SiO2=y_SiO2_Au3nm_SiO2{12,1};

%x_Fermi_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{1,1};
%x_Elvl10_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{2,1};
x_Elvl9_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{3,1};
x_Elvl8_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{4,1};
x_Elvl7_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{5,1};
%x_Elvl6_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{6,1};
%x_Elvl5_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{7,1};
%x_Elvl4_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{8,1};
%x_Elvl3_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{9,1};
%x_Elvl2_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{10,1};
%x_Elvl1_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{11,1};
%x_Well_SiO2_Au3nm_SiO2=x_SiO2_Au3nm_SiO2{12,1};


%% SiO2 Au Al2O3 Wells

h_SiO2_Au3nm_HfO2=openfig('/Users/mainuser/Desktop/Au SHG Paper Figures/SiO2_Au3nm_Al2O3_with_fermi_energy.fig');
h_SiO2_Au3nm_HfO2=findobj(gca,'Type','Line');
x_SiO2_Au3nm_Al2O3=get(h_SiO2_Au3nm_HfO2,'Xdata');
y_SiO2_Au3nm_Al2O3=get(h_SiO2_Au3nm_HfO2,'Ydata');

%y_Fermi_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{1,1};
y_Elvl9_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{2,1};
y_Elvl8_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{3,1};
y_Elvl7_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{4,1};
%y_Elvl6_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{5,1};
%y_Elvl5_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{6,1};
%y_Elvl4_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{7,1};
%y_Elvl3_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{8,1};
%y_Elvl2_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{9,1};
%y_Elvl1_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{10,1};
%y_Well_SiO2_Au3nm_Al2O3=y_SiO2_Au3nm_Al2O3{11,1};

%x_Fermi_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{1,1};
x_Elvl9_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{2,1};
x_Elvl8_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{3,1};
x_Elvl7_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{4,1};
%x_Elvl6_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{5,1};
%x_Elvl5_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{6,1};
%x_Elvl4_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{7,1};
%x_Elvl3_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{8,1};
%x_Elvl2_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{9,1};
%x_Elvl1_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{10,1};
%x_Well_SiO2_Au3nm_Al2O3=x_SiO2_Au3nm_Al2O3{11,1};

%% SiO2 Au HfO2 Wells

h_SiO2_Au3nm_HfO2=openfig('/Users/mainuser/Desktop/Au SHG Paper Figures/SiO2_Au3nm_HfO2_with_fermi_energy.fig');
h_SiO2_Au3nm_HfO2=findobj(gca,'Type','Line');
x_SiO2_Au3nm_HfO2=get(h_SiO2_Au3nm_HfO2,'Xdata');
y_SiO2_Au3nm_HfO2=get(h_SiO2_Au3nm_HfO2,'Ydata');

%y_Fermi_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{1,1};
y_Elvl9_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{2,1};
y_Elvl8_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{3,1};
y_Elvl7_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{4,1};
%y_Elvl6_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{5,1};
%y_Elvl5_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{6,1};
%y_Elvl4_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{7,1};
%y_Elvl3_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{8,1};
%y_Elvl2_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{9,1};
%y_Elvl1_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{10,1};
%y_Well_SiO2_Au3nm_HfO2=y_SiO2_Au3nm_HfO2{11,1};

%x_Fermi_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{1,1};
x_Elvl9_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{2,1};
x_Elvl8_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{3,1};
x_Elvl7_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{4,1};
%x_Elvl6_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{5,1};
%x_Elvl5_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{6,1};
%x_Elvl4_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{7,1};
%x_Elvl3_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{8,1};
%x_Elvl2_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{9,1};
%x_Elvl1_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{10,1};
%x_Well_SiO2_Au3nm_HfO2=x_SiO2_Au3nm_HfO2{11,1};


%%
%EOF