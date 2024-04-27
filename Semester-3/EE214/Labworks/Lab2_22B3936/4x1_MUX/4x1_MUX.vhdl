library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity MUX_4x1 is
port ( In4, In3, In2, In1, S2, S1: in std_logic; Y0: out std_logic);
end entity MUX_4x1;

architecture Struct of MUX_4x1 is
signal MUX_out1, MUX_out2 : std_logic;

component MUX_2x1 is
     port(in1,in0,S: in std_logic; Y0: out std_logic);
   end component;

begin
  --component instances
  MUX1 : MUX_2x1 port map (in1 => In2, in0 => In1, S => S1, Y0 => MUX_out1);
  MUX2 : MUX_2x1 port map (in1 => In4, in0 => In3, S => S1, Y0 => MUX_out2);
  MUX3 : MUX_2x1 port map (in1 => MUX_out2, in0 => MUX_out1, S => S2, Y0 => Y0);
  end Struct;