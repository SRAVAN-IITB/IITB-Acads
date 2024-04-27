library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity MUX_4x1_4bit is
port ( A, B, C, D: in std_logic_vector(3 downto 0); S: in std_logic_vector(1 downto 0); Y0: out std_logic_vector(3 downto 0));
end entity MUX_4x1_4bit;

architecture Struct of MUX_4x1_4bit is

--component MUX_2x1 is
--     port(in1,in0,S: in std_logic; Y0: out std_logic);
--   end component;
	
component MUX_4x1 is
port ( In4, In3, In2, In1, S2, S1: in std_logic; Y0: out std_logic);
end component;

begin
  --component instances
  MUX1 : MUX_4x1 port map (In1 => A(0), In2 => B(0), In3 => C(0), In4 => D(0), S2 => S(1), S1 => S(0), Y0 => Y0(0));
  MUX2 : MUX_4x1 port map (In1 => A(1), In2 => B(1), In3 => C(1), In4 => D(1), S2 => S(1), S1 => S(0), Y0 => Y0(1));
  MUX3 : MUX_4x1 port map (In1 => A(2), In2 => B(2), In3 => C(2), In4=> D(2), S2 => S(1), S1 => S(0), Y0 => Y0(2));
  MUX4 : MUX_4x1 port map (In1 => A(3), In2 => B(3), In3 => C(3), In4 => D(3), S2 => S(1), S1 => S(0), Y0 => Y0(3));
  end Struct;