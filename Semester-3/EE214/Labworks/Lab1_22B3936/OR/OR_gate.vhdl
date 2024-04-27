library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity OR_gate is
port (A, B: in std_logic; Y: out std_logic);
end entity OR_gate;

architecture Struct of OR_gate is
signal mid : std_logic;
begin
  --component instances
  NOR1 : NOR_2 port map (A => A, B => B, Y => mid);
  -- FINAL GATE
  NOR2 : NOR_2 port map (A=>mid, B=>mid, Y=>Y);
  end Struct;