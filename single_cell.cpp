#include "config.hpp"
#include "all_particles.hpp"
#include "Cell.hpp"
#include "utils/for_each_pair.hpp"
#include "velocity_verlet_inline.hpp"
#include "timer.hpp"

#include <iostream>
#include <string>



template <typename Particle>
Cell<Particle> setup_cell(int n) {
Cell<Particle> cell;
for (int i=0;i<n;i++) cell.particles().insert(Particle{});
return cell; 
};

template <typename ParticleIterable>
void all_with_all(ParticleIterable& particles) {
  auto fake_force_kernel =[](auto &p1, auto& p2) {
    auto d =p1.pos()-p2.pos();
    p1.force()+=d;
    p2.force() -=d;
  };
  
  Utils::for_each_pair(particles.begin(), particles.end(), fake_force_kernel);
};


template <typename ParticleIterable> 
void single_pass(ParticleIterable& particles) {
  velocity_verlet_step_1(particles, 0.01);
}


template <typename Particle>
void time(int n) {
  auto cell = setup_cell<Particle>(n);
  auto tick = Timer::now();
  single_pass(cell.particles());
  auto single_pass_time = Duration(Timer::now() -tick);
  
  tick = Timer::now();
  all_with_all(cell.particles());
  auto all_with_all_time = Duration(Timer::now() -tick);
//  std::cerr << (cell.particles().begin())->pos() <<std::endl;
  std::cout <<n <<" " <<to_ms(single_pass_time) << " "<<to_ms(all_with_all_time)<<std::endl;
}


template <typename Particle>
void run_timing(std::string particle_name) {
  std::cout << particle_name <<" Size: "<<sizeof(Particle)<< std::endl;
  for (int n: {5000, 10000, 15000, 20000,25000,30000,35000,40000, 45000, 50000, 55000, 60000, 65000, 70000, 75000, 80000}) {
    time<Particle>(n);
  }
}
int main(int argc, char** argv) {
  run_timing<CurrentEsParticle>("Current Espresso Particle");
  run_timing<MinimalFlatParticle<456>>("Re-ordered properties");
  run_timing<MinimalFlatParticle<0>>("Minimal flat particle");
  run_timing<SoABackedParticle>("Struct-of-arrays backed particle");
  return 0;
}
