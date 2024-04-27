library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity AND_gate is
port (A, B: in std_logic; Y: out std_logic);
end entity AND_gate;

architecture Struct of AND_gate is
signal A_BAR,B_BAR : std_logic;
begin
  --component instances
  NOR1 : NOR_2 port map (A => A, B => A, Y => A_BAR);
  NOR2 : NOR_2 port map (A => B, B => B, Y => B_BAR);
  -- FINAL GATE
  NOR3 : NOR_2 port map (A=>A_BAR, B=>B_BAR, Y=>Y);
  end Struct;