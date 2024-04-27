library ieee;
use ieee.std_logic_1164.all;

library work;

entity SequenceGenerator is
	port (reset,clock: in std_logic; y:out std_logic);
end entity SequenceGenerator;

architecture struct of SequenceGenerator is
	signal D :std_logic_vector(2 downto 0);
	signal Q:std_logic_vector(2 downto 0);
	
	begin
		
		d_flip_flops: process (clock,reset)
			begin
			if(reset='1')then -- set implies flip flip output logic high
			Q(0) <= '1';
			Q(1) <= '0';
			Q(2) <= '0';
			elsif (clock'event and (clock='1')) then
			Q(0) <= D(0);
			Q(1)<= D(1);
			Q(2) <= D(2);
			end if ;
			
			if (Q(2)='0' and Q(1)='0' and Q(0)='1') then
				D<="011";
			elsif (Q(2)='0' and Q(1)='1' and Q(0)='1') then
				D<="101";
			elsif (Q(2)='1' and Q(1)='0' and Q(0)='1') then
				D<="000";
			elsif (Q(2)='0' and Q(1)='0' and Q(0)='0') then
				D<="001";
			else
				D<="001";
			end if;
		end process d_flip_flops;
		
--		D(2)<=Q(1) and Q(0) and (not Q(2));
--		D(1)<= (not Q(2) ) and (not Q(1) ) and Q(0);
--		D(0)<= not ( Q(2) and (not Q(1)) and Q(0));
		
		y<=Q(0);
end struct;