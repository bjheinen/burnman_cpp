#include <iostream>
#include <string>
#include "burnman/core/mineral.hpp"
#include "burnman/utils/eos.hpp"

const char* to_string(EOSType eos_type) {
    switch (eos_type) {
        case EOSType::BM2: return "BM2";
        case EOSType::BM3: return "BM3";
        case EOSType::MGD3: return "MGD3";
        default: return "?";
    }
}

int main()
{
  
  std::cout << "Making a mineral class!\n";
  Mineral fp;
  std::cout << "Get name (should class name): " << fp.get_name() << "\n";
  fp.set_name("Fe-Periclase");
  std::cout << "Get name (should be real): " << fp.get_name() << "\n";

  std::cout << "Setting parameters\n";
  fp.params.V_0 = 22.9e-6;
  fp.params.K_0 = 157.5e9;
  fp.params.Kprime_0 = 3.92;
  fp.params.molar_mass = 0.04567;
  fp.params.napfu = 2;
  fp.params.debye_0 = 587;
  fp.params.grueneisen_0 = 1.46;
  fp.params.q_0 = 1.2;
  std::cout << "Parameters set\n";
  std::cout << "Trying to grab molar mass\n";
  std::cout << "Molar mass = " << fp.get_molar_mass() << "\n";
  std::cout << "Setting EOS\n";
  fp.set_method(EOSType::MGD3);
  std::cout << "EOS type = " << to_string(*fp.params.equation_of_state) << "\n";
  std::cout << "Setting P & T to 2 GPa and 873 K\n";
  fp.set_state(2.e9, 873.0);
  std::cout << "Querying some properties\n";

  std::cout << "Pressure = " << fp.get_pressure() << "\n";
  std::cout << "Temperature = " << fp.get_temperature() << "\n";
  std::cout << "Internal energy = " << fp.get_molar_internal_energy() << "\n";
  std::cout << "Gibbs = " << fp.get_molar_gibbs() << "\n";
  std::cout << "Helmholtz = " << fp.get_molar_helmholtz() << "\n";
  std::cout << "Mol Mass = " << fp.get_molar_mass() << "\n";
  std::cout << "Volume = " << fp.get_molar_volume() << "\n";
  std::cout << "Unmod vol = " << fp.get_molar_volume_unmodified() << "\n";
  std::cout << "Density = " << fp.get_density() << "\n";
  std::cout << "Entropy = " << fp.get_molar_entropy() << "\n";
  std::cout << "Enthalpy = " << fp.get_molar_enthalpy() << "\n";
  std::cout << "K_T = " << fp.get_isothermal_bulk_modulus_reuss() << "\n";
  std::cout << "K_S = " << fp.get_isentropic_bulk_modulus_reuss() << "\n";
  std::cout << "1/K_T = " << fp.get_isothermal_compressibility_reuss() << "\n";
  std::cout << "1/K_S = " << fp.get_isentropic_compressibility_reuss() << "\n";
  std::cout << "G = " << fp.get_shear_modulus() << "\n";
  std::cout << "Vp = " << fp.get_p_wave_velocity() << "\n";
  std::cout << "Vphi = " << fp.get_bulk_sound_velocity() << "\n";
  std::cout << "Vs = " << fp.get_shear_wave_velocity() << "\n";
  std::cout << "gamma = " << fp.get_grueneisen_parameter() << "\n";
  std::cout << "alpha = " << fp.get_thermal_expansivity() << "\n";
  std::cout << "Cv = " << fp.get_molar_heat_capacity_v() << "\n";
  std::cout << "Cp = " << fp.get_molar_heat_capacity_p() << "\n";
  std::cout << "therm grad = " << fp.get_isentropic_thermal_gradient() << "\n";

  std::cout << "Yes?" << std::endl;



}