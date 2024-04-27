library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Flipflops.all;
use work.Gates.all;

entity Sequence_generator_stru_dataflow is
port (reset,clock: in std_logic;
y:out std_logic);
end entity Sequence_generator_stru_dataflow;
architecture struct of Sequence_generator_stru_dataflow is
signal D:std_logic_vector(2 downto 0);
signal Q:std_logic_vector(2 downto 0);

begin
    -- equations for D(2), D(1), and D(0)
D(2) <= (Q(2) and (not Q(1)) AND Q(0)) or ((not Q(0)) and  Q(1)AND (NOT Q(2)));
D(1) <= (Q(0) and (not Q(1))) or ((not Q(2)) AND (not Q(1)));
D(0) <= Q(2) or (Q(0) xor Q(1));

DFF2 : dff_reset port map (D => D(2), clock => clock, reset => reset, Q => Q(2));
DFF1 : dff_reset port map (D => D(1), clock => clock, reset => reset, Q => Q(1));
DFF0 : dff_set port map (D => D(0), clock => clock, set => reset, Q => Q(0));

y <= Q(0);

-- write the equations in dataflow e.g z=a+bc written as z <= a or (b and c)
-- for DFF inputs which was derived and also for y.
-- Instantiate components dff_reset
-- and dff_set appropriately using port map statements.
end struct;