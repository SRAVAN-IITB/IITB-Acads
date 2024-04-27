library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity sequence is
	port(reset,clock:in std_logic;
			outp: out std_logic);
end sequence;

architecture bhv of sequence is
---------------Define state type here-----------------------------
	type state is (s0,s1,s2,s3); -- Fill other states here
---------------Define signals of state type-----------------------
	signal s: state:=s0;
	begin
		clock_proc:process(clock,reset)
		begin
			if(reset = '1') then
				s <= s0;
			elsif(clock'event and clock = '1') then
				 case s is
					  when s0 =>
							s <= s1;
					  when s1 =>
							s <= s2;
					  when s2 =>
							s <= s3;
					  when s3 =>
							s <= s0;
				 end case;
			end if;
		end process;

		outpproc_proc:process(clock, s)
		begin
            case s is
                when s0 =>
                    outp <= '1';
                when s1 =>
                    outp <= '1';
                when s2 =>
                    outp <= '1';
                when s3 =>
                    outp <= '0';
            end case;
        end process;
end bhv;