library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity Servers is
	port(clock:in std_logic;
        input: in std_logic_vector(3 downto 0);
		outp: out std_logic_vector(4 downto 0));
end Servers;

architecture bhv of Servers is
---------------Define state type here-----------------------------
	type state is (idle, Ad, Pr, Re, St, PAW, RAW, RPW, SRW, SPW, SAW); -- Fill other states here
---------------Define signals of state type-----------------------
	signal s: state:=idle;
    signal n: integer := 1;
	begin
		clock_proc:process(clock, s)
		begin
        if(clock'event and clock = '1') then
                case s is
                    when idle =>
                        if(input(3)= '1') then
                            s <= Ad;
									 n <= 1;
                        elsif(input(2) = '1') then
                            s <= Pr;
									 n <= 1;
                        elsif(input(1) = '1') then
                            s <= Re;
									 n <= 1;
								elsif(input(0) = '1') then
									 s <= St;
									 n <= 1;
                        else
                            s <= idle;
									 n <= 1;
                        end if;

                    when Ad =>
                        if(input(3) = '1') then
                            s <= Ad;
									 n <= 1;
                        else
                            s <= idle;
									 n <= 1;
                        end if;
                    
                    when Pr =>
                        if(input(3) = '1') then
                            s <= PAW;
                        elsif(input(2) = '1') then
                            s <= Pr;
                        else
                            s <= idle;
                        end if;
								n <= 1;
                    
                    when Re =>
                        if(input(3) = '1') then
                            s <= RAW;
                        elsif(input(2) = '1') then
                            s <= RPW;
                        elsif(input(1) = '1') then
                            s <= Re;
                        else
                            s <= idle;
                        end if;
								n <= 1;

                    when St =>
                        if(input(3) = '1') then
                            s <= SAW;
                        elsif(input(2) = '1') then
                            s <= SPW;
                        elsif(input(1) = '1') then
                            s <= SRW;
								elsif(input(0) = '1') then
									 s <= St;
                        else
                            s <= idle;
                        end if;
								n <= 1;

                    when PAW =>
                        if(input(3) = '1') then
                            if(n = 150000000) then
                                n <= 1;
                                s <= Ad;
                            else
                                n <= n + 1;
                                s <= PAW;
                            end if;
                        elsif(input(2) = '1') then
                            s <= Pr;
                            n <= 1;
                        else
                            s <= idle;
                            n <= 1;
                        end if;

                    when RAW =>
                        if(input(3) = '1') then
                            if(n = 150000000) then
                                n <= 1;
                                s <= Ad;
                            else
                                n <= n + 1;
                                s <= RAW;
                            end if;
                        elsif(input(2) = '1') then
                            s <= RPW;
                            n <= n + 1;
                        elsif(input(1) = '1') then
                            s <= Re;
                            n <= 1;
                        else
                            s <= idle;
                            n <= 1;
                        end if;
                    
                    when RPW =>
                        if(input(3) = '1') then
                            s <= RAW;
                            n <= n + 1;
                        elsif(input(2) = '1') then
                            if(n = 150000000) then
                                n <= 1;
                                s <= Pr;
                            else
                                n <= n + 1;
                                s <= RPW;
                            end if;
                        elsif(input(1) = '1') then
                            s <= Re;
                            n <= 1;
                        else
                            s <= idle;
                            n <= 1;
                        end if;    
                    
                    when SAW =>
                        if(input(3) = '1') then
                            if(n = 150000000) then
                                n <= 1;
                                s <= Ad;
                            else
                                n <= n + 1;
                                s <= SAW;
                            end if;
                        elsif(input(2) = '1') then
                            s <= SPW;
                            n <= n + 1;
                        elsif(input(1) = '1') then
                            s <= SRW;
                            n <= n + 1;
                        elsif(input(0) = '1') then
                            s <= St;
                            n <= 1;
                        else
                            s <= idle;
                            n <= 1;
                        end if;
                    
                    when SPW =>
                        if(input(3) = '1') then
                            s <= SAW;
                            n <= n + 1;
                        elsif(input(2) = '1') then
                            if(n = 150000000) then
                                n <= 1;
                                s <= Pr;
                            else
                                n <= n + 1;
                                s <= SPW;
                            end if;
                        elsif(input(1) = '1') then
                            s <= SRW;
                            n <= n + 1;
                        elsif(input(0) = '1') then
                            s <= St;
                            n <= 1;
                        else
                            s <= idle;
                            n <= 1;
                        end if;
                    
                    when SRW =>
                        if(input(3) = '1') then
                            s <= SAW;
                            n <= n + 1;
                        elsif(input(2) = '1') then
                            s <= SPW;
                            n <= n + 1;
                        elsif(input(1) = '1') then
                            if(n = 150000000) then
                                n <= 1;
                                s <= Re;
                            else
                                n <= n + 1;
                                s <= SRW;
                            end if;
                        elsif(input(0) = '1') then
                            s <= St;
                            n <= 1;
                        else
                            s <= idle;
                            n <= 1;
                        end if;
							when others =>
								s <= idle;
								n <= 1;
                end case;
        end if;
		end process;

		outpproc_proc:process(clock, s)
		begin
            case s is
					 when idle =>
						  outp <= "00000";
                when Ad =>
                    outp <= "01111";
                when Pr =>
                    outp <= "00111";
                when Re=>
                    outp <= "00011";
                when St=>
                    outp <= "00001";
                when others=>
                    outp <= "11111";
            end case;
        end process;
end bhv;