library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity prime_det is
port (in4, in3, in2, in1, in0: in std_logic; Y0: out std_logic);
end entity prime_det;

architecture Struct of prime_det is
signal x0orx3, nx2, nx2ANDx1, x2103, nx1, nx4, nx3, and210, orringfor_nx4, x_3210, 
			orred_nx4, xnor32, anded3 : std_logic;

begin
  --component instances
  x1inv : INVERTER port map (A=> in1, Y => nx1); 
  x2inv : INVERTER port map (A=> in2, Y => nx2);  
  x3inv : INVERTER port map (A=> in3, Y => nx3);
  x4inv : INVERTER port map (A=> in4, Y => nx4);
  x03or : OR_2 port map (A => in3, B => in4, Y => x0orx3);
  XNOR1 : XNOR_2 port map (A => in3, B => in2, Y => xnor32);
  x2inv_x1AND : AND_2 port map (A=> nx2, B => in1, Y => nx2ANDx1);
  AND2103 : AND_2 port map (A=> nx2ANDx1, B => x0orx3, Y => x2103);
  AND2n10 : AND_3 port map (A=> in2, B => in0, C => nx1, Y => and210);
  orringforX4bar : OR_2 port map (A=> and210, B => x2103, Y => orringfor_nx4);
  OR3210 : OR_4 port map (A=> nx3, B => in2, C => in1, D => in0, Y => x_3210);
  orringX4bar : OR_2 port map (A=> nx4, B => orringfor_nx4, Y => orred_nx4);
  anding3 : AND_3 port map (A=> in4, B => in0, C => xnor32, Y => anded3);
  finalORRING: OR_3 port map (A=> anded3, B => orred_nx4, C => x_3210, Y => Y0); 

  end Struct;