library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity MUX_2x1 is
port (in1, in0, S: in std_logic; Y0: out std_logic);
end entity MUX_2x1;

architecture Struct of MUX_2x1 is
signal AND_out1, AND_out2, S_bar : std_logic;

begin
  --component instances
  NOT1 : INVERTER port map (A => S, Y => S_bar);
  AND1 : AND_2 port map (A => in0, B => S_bar, Y => AND_out1);
  AND2 : AND_2 port map (A=> in1, B => S, Y => AND_out2);
  OR1 : OR_2 port map (A=> AND_out1, B => AND_out2, Y => Y0);
  end Struct;