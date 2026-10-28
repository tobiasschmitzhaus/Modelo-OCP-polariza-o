import numpy as np

class MotorTermodinamico:
    def __init__(self):
        pass 

    def calcular_pH_alcalinidade(self, AT_m3, DIC_m3, I_m3, Fe2_total_m3):
        AT = AT_m3 / 1000.0
        DIC = np.maximum(DIC_m3 / 1000.0, 1e-16)
        Fe2_t = np.maximum(Fe2_total_m3 / 1000.0, 1e-16)
        I = np.clip(I_m3 / 1000.0, 1e-5, 0.5)
        
        A = 0.509
        sqrt_I = np.sqrt(I)
        log_gamma_1 = -A * 1.0**2 * (sqrt_I / (1.0 + sqrt_I) - 0.3 * I)
        log_gamma_2 = -A * 2.0**2 * (sqrt_I / (1.0 + sqrt_I) - 0.3 * I)
        gamma_1 = 10**log_gamma_1
        gamma_2 = 10**log_gamma_2

        Kw = 1e-14 / (gamma_1**2)
        Ka1 = (10**-6.35) / (gamma_1**2)
        Ka2 = (10**-10.33) / gamma_2
        Kh1 = (10**-9.5) * gamma_2 / (gamma_1**2)
        Ksp = 8.0e-16

        # Busca no domínio de pH (0 a 14) para eliminar o teto numérico da bisseção linear
        pH_min = np.full_like(AT, 1.0)
        pH_max = np.full_like(AT, 14.0)

        for _ in range(35): 
            pH_mid = (pH_min + pH_max) / 2.0
            H_mid = 10**(-pH_mid)
            
            # Sem limite numérico arbitrário. O OH- livre dita o balanço no líquido.
            # A precipitação já é cuidada pelo r_precip_bulk na cinética.
            OH_mid = Kw / H_mid
            
            denom_c = H_mid**2 + Ka1 * H_mid + Ka1 * Ka2
            HCO3 = DIC * (Ka1 * H_mid) / denom_c
            CO3  = DIC * (Ka1 * Ka2) / denom_c
            
            # Conservação EXATA da massa do Ferro (Evita que FeOH > Fe2_total)
            Fe_livre = Fe2_t / (1.0 + Kh1 / H_mid)
            FeOH = Fe_livre * Kh1 / H_mid
            
            F_mid = OH_mid - H_mid + HCO3 + 2.0 * CO3 + FeOH - AT
            
            pH_min = np.where(F_mid < 0, pH_mid, pH_min)
            pH_max = np.where(F_mid >= 0, pH_mid, pH_max)

        pH_eq = (pH_min + pH_max) / 2.0
        H_eq = 10**(-pH_eq)
        OH_eq = Kw / H_eq
        return H_eq * 1000.0, OH_eq * 1000.0

    def calcular_SI_FeOH2(self, Fe2_m3, H_m3, OH_m3):
        Fe2 = np.maximum(Fe2_m3 / 1000.0, 1e-16)
        OH = np.maximum(OH_m3 / 1000.0, 1e-16)
        
        I = np.clip(3.0 * Fe2 + 0.1, 1e-5, 3.0)
        A = 0.509
        sqrt_I = np.sqrt(I)
        log_gamma_1 = -A * 1.0**2 * (sqrt_I / (1.0 + sqrt_I) - 0.3 * I)
        log_gamma_2 = -A * 2.0**2 * (sqrt_I / (1.0 + sqrt_I) - 0.3 * I)
        
        gamma_1 = 10**log_gamma_1 
        gamma_2 = 10**log_gamma_2 
        Ksp = 8.0e-16
        
        log10_Q = np.log10(gamma_2 * Fe2) + 2.0 * np.log10(gamma_1 * OH)
        SI = log10_Q - np.log10(Ksp)
        
        return np.clip(SI, -20.0, 20.0)