library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity bitshifter is
port (a: in std_logic_vector(7 downto 0); 
		b: in std_logic_vector(2 downto 0); 
		L: in std_logic; 
		s: out std_logic_vector(7 downto 0);
end entity bitshifter;

architecture Struct of bitshifter is
signal s1, s2, s3, s4: std_logic_vector(7 downto 0);

component MUX_2x1 is
port (in1, in0, S: in std_logic; Y0: out std_logic);
end component MUX_2x1;

begin
  --component instances
  bitrev1 : INVERTER port map (A => S, Y => S_bar);
  AND1 : AND_2 port map (A => in0, B => S_bar, Y => AND_out1);
  AND2 : AND_2 port map (A=> in1, B => S, Y => AND_out2);
  OR1 : OR_2 port map (A=> AND_out1, B => AND_out2, Y => Y0);
  end Struct;