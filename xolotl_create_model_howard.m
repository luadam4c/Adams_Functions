function xolotlObject = xolotl_create_model_howard
%% Creates a xolotl model based on Howard et al
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 2019-08-16 Adapted from Alecs_Functions/xolotl-model-zoo/model_howard.m
% 2019-08-16 Renamed compartment as 'soma'

xolotlObject = xolotl('temperature', 36, 'temperature_ref', 22);

xolotlObject.add('compartment', 'soma', 'Cm', 10, 'A', 1e-3, 'vol', 1e-6, 'Ca_out', 2e3);

xolotlObject.soma.add('buchholtz/CalciumMech', 'Ca_in', 0.2, 'tau_Ca', 250, 'phi', 1);

xolotlObject.soma.add('soplata/CaN', 'gbar', 0.25);
xolotlObject.soma.add('destexhe/CaT', 'gbar', 1);
xolotlObject.soma.add('soplata/Kd', 'gbar', 200);
xolotlObject.soma.add('soplata/NaV', 'gbar', 2000);
xolotlObject.soma.add('Leak', 'gbar', 0.1);

xolotlObject.I_ext = 0.07;
