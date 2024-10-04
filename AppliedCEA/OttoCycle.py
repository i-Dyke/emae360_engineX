import copy as c 
import CEA_Wrap as cea

class ThermoState: 
    def __init__(self, pressure_kPa:float, temperature_K:float, density_kgm3:float, volume_m3:float):
        self.pressure_kPa = pressure_kPa
        self.temperature_K = temperature_K
        self.density_kgm3 = density_kgm3
        self.volume_m3 = volume_m3
            
    def set_pressure(self, p_kPa:float):
        self.pressure_kPa =  p_kPa
    
    def set_temperature(self, T_K:float):
        self.temperature_K = T_K
    
    def set_density(self, d_kgm3:float):
        self.density_kgm3 = d_kgm3
    
    def set_volume(self, v_m3:float):
        self.volume_m3 = v_m3
    
    def get_pressure(self): 
        return self.pressure_kPa
    def get_temperature(self): 
        return self.temperature_K
    def get_density(self): 
        return self.density_kgm3
    def get_volume(self): 
        return self.volume_m3

    def curr_state(self):
        vals = [self.pressure_kPa, self.temperature_K, self.density_kgm3, self.volume_m3]
        vals = [str(i) for i in vals]
        print('\t'.join(["","p kPa", vals[0], '\n', "T K", vals[1], '\n', "d kg/m3", vals[2], '\n', "v m3", vals[3]]))

class Otto:
    def __init__(self, p_kPa:float, T_K:float, d_kgm3:float, v_m3:float, CR:float, gamma:float, comp_eff:float, exp_eff:float):
        self.state = ThermoState(p_kPa, T_K, d_kgm3, v_m3)
        self.CR= CR
        self.gamma= gamma
        self.comp_eff= comp_eff
        self.exp_eff= exp_eff

    # State1 
    def intake(self):
        self.state.set_volume(self.state.get_volume() * self.CR)

    #State2
    def compression(self):
        self.state.set_volume(v_m3= self.state.get_volume() / self.CR)
        self.state.set_density(self.state.get_density() * self.CR)
        self.state.set_pressure(self.state.get_pressure() * self.CR**(self.gamma))
        self.state.set_temperature(T_K= self.state.get_temperature() * ( 1 + (self.CR**(self.gamma - 1) -1)/self.comp_eff ))

    #State3
    def combustion_ideal(self, IGCon:float, QComb:float, Cp:float, comb_eff:float, o_f:float):
        self.state.set_temperature(T_K= self.state.get_temperature() + comb_eff*(QComb / (Cp * o_f)))
        # Note: now using post-combustion temp
        self.state.set_pressure(p_kPa= self.state.get_density() * self.state.get_temperature() * IGCon)

    #State4
    def powStroke(self):
        self.state.set_density(d_kgm3= self.state.get_density() / self.CR)
        self.state.set_volume(v_m3= self.state.get_volume() * self.CR)
        self.state.set_temperature(T_K= self.state.get_temperature() * (1 - self.exp_eff*(1 - self.CR**(1-self.gamma))))
        self.state.set_pressure(p_kPa= self.state.get_pressure() * (1/self.CR)**(self.gamma))

    #State5
    def exhaust(self, orig:ThermoState):
        self.state.set_density(orig.get_density())
        self.state.set_volume(orig.get_volume())
        self.state.set_temperature(orig.get_temperature())
        self.state.set_pressure(orig.get_pressure())  
    #State6

    def get_state(self):
        return self.state

def cyc_ideal(Otto:Otto, IGCon:float, QComb:float, Cp:float, comb_eff:float, o_f:float):

    state= Otto.get_state()
    # Pre-Intake
    s1= c.deepcopy(state)

    Otto.intake()
    # Pre-Compression
    s2= c.deepcopy(state)
        
    Otto.compression()
    # Pre-Combustion
    s3= c.deepcopy(state)
        
    Otto.combustion_ideal(IGCon, QComb, Cp, comb_eff, o_f)
    # Post-Combustion
    s4 = c.deepcopy(state)
        
    Otto.powStroke()
    # Pre-Exhaust
    s5 = c.deepcopy(state)
        
    Otto.exhaust(s1)
    # Post-Exhaust/Pre-Intake
    s6 = c.deepcopy(state)

    return [s1,s2,s3,s4,s5,s6]


def cyc_real(Otto:Otto, File:str, Mat:cea.Material, massf:bool, o_f:float, reac_warming:bool=False, reac_temp:float=298.15 
, report_comb:bool=False):
    
    state= Otto.get_state()
    # Pre-Intake
    s1= c.deepcopy(state)

    Otto.intake()
    # Pre-Compression
    s2= c.deepcopy(state)
        
    Otto.compression()
    # Pre-Combustion
    s3= c.deepcopy(state)

    for i in Mat:
        if reac_warming:
            i.set_temp(state.get_temperature())
        else:
            i.set_temp(reac_temp)
    #print(' '.join([str(i.temp) for i in Mat]))

    combProb = cea.UVProblem(
        filename = File,
        materials = Mat,
        density = state.get_density(),
        density_units = "kg",
        massf = massf,
        o_f = o_f
    )
    realComb = combProb.run()
    # Inputting Real Combustion into the State Directly
    state.set_temperature(realComb.t)
        # *100 to turn from bar to kPa
    state.set_pressure(realComb.p * 100)

    # Post-Combustion
    s4 = c.deepcopy(state)
        
    Otto.powStroke()
    # Pre-Exhaust
    s5 = c.deepcopy(state)
        
    Otto.exhaust(s1)
    # Post-Exhaust/Pre-Intake
    s6 = c.deepcopy(state)

    if report_comb:
        return [s1,s2,s3,s4,s5,s6,realComb]
    else:
        return [s1,s2,s3,s4,s5,s6]

def cyc_work(Cycle_States:ThermoState, CR:float, Cv_kJ:float, CylCt:int, Disp_m3:float):
    if len(Cycle_States) < 6:
        raise Exception("Require a cycle list input of 6+ states")
    
    cs = Cycle_States
    return (Cv_kJ * (cs[3].get_temperature() -cs[2].get_temperature()) - Cv_kJ * (cs[4].get_temperature() -cs[1].get_temperature())) * (CylCt * Disp_m3)*(CR/(CR-1)) * cs[1].get_density() 

    