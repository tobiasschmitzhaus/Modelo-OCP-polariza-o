import numpy as np

class MotorCinetico:
    def __init__(self):
        self.F = 96485.33
        self.R = 8.314
        self.T = 298.15
        
        self.C_dl = 0.2  # F/m2
        
        # =================================================================
        # POTENCIAIS DE EQUILÍBRIO PADRÃO (SHE @ pH 0 ou estados base)
        # =================================================================
        self.E_eq_Fe_padrao  = -0.44   # Fe/Fe2+
        self.E_eq_O2_puro    = 1.229   # O2 + 4H+ + 4e- -> 2H2O
        self.E_eq_H2O_puro   = 0.0     # 2H+ + 2e- -> H2
        self.E_eq_NO3_padrao = 0.01    # NO3- + H2O + 2e- -> NO2- + 2OH-
        
        # =================================================================
        # CORRENTES DE TROCA BASE (A/m²) E CONSTANTES CINÉTICAS
        # Recalibrado para eletrólito 0.1 M NaNO3 (base empírica)
        # =================================================================
        self.i_0_Fe     = 1.0e-2   # Escalonado para o estado ativo massivo (icorr ~ 0.7 A/m2)
        self.i_0_O2_ref = 5.0e-4   # Ajustado para balancear a nova cinética
        self.i_0_H2O    = 1.0e-9   
        self.i_0_H      = 1.0e-6   
        self.i_0_NO3    = 1.0e-4   # Base para a redução direta do nitrato
        self.k_crust    = 5.0e-3   # Constante de redução da goethita para magnetita
        
        self.a_OH_ref = 1e-4     
        self.ordem_OH = 1.0      
        
        # =================================================================
        # COEFICIENTES DE TRANSFERÊNCIA (b = 2.303*RT / beta*F)
        # Derivados de b_a = 262 mV/dec e b_c = 98 mV/dec
        # =================================================================
        self.beta_a_Fe  = 0.226  # Corrigido da extrapolação anódica
        self.beta_c_Fe  = 0.5
        self.beta_c_O2  = 0.604  # Corrigido da extrapolação catódica
        self.beta_c_H2O = 0.5
        self.beta_c_NO3 = 0.5
        self.beta_c_crust = 0.5
        
        # =================================================================
        # DIFUSIVIDADE E TRANSPORTE
        # =================================================================
        self.D_Na = 1.33e-9
        self.z_Na = 1.0
        self.D_NO3 = 1.90e-9
        self.z_NO3 = -1.0

        # =================================================================
        # SATURAÇÃO DE OXIGÊNIO (Lei de Henry)
        # =================================================================
        self.K_H_O2 = 0.0013  # Constante de Henry para O2 em mol/(L atm) a 298K
        self.P_O2 = 0.21      # Pressão parcial do oxigênio na atmosfera padrão
        self.C_O2_sat = (self.K_H_O2 * self.P_O2) * 1000.0  # Conversão de mol/L para mol/m3