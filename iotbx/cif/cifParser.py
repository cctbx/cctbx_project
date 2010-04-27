# $ANTLR 3.1.2 cif.g 2010-04-27 15:19:51

import sys
from antlr3 import *
from antlr3.compat import set, frozenset

import iotbx.cif
from cctbx.array_family import flex



# for convenience in actions
HIDDEN = BaseRecognizer.HIDDEN

# token types
DOUBLE_QUOTED_STRING=29
CHAR_STRING=13
EXPONENT=12
NON_BLANK_CHAR=27
LOOP_HEADER=9
SEMI_COLON_TEXT_FIELD=14
SINGLE_QUOTED_STRING=28
DOUBLE_QUOTE=16
GLOBAL_=24
ORDINARY_CHAR=18
WHITESPACE=5
VERSION=26
EOF=-1
TAG=8
SINGLE_QUOTE=17
T__31=31
T__32=32
STOP_=25
EOL=15
T__33=33
NON_BLANK_CHAR_=19
T__34=34
T__35=35
COMMENTS=4
T__36=36
SAVE_FRAME_HEADING=6
ANY_PRINT_CHAR=21
TEXT_LEAD_CHAR=20
SAVE_=7
LOOP_=23
DIGIT=11
UNQUOTED_STRING=30
DATA_=22
DATA_BLOCK_HEADING=10

# token names
tokenNames = [
    "<invalid>", "<EOR>", "<DOWN>", "<UP>",
    "COMMENTS", "WHITESPACE", "SAVE_FRAME_HEADING", "SAVE_", "TAG", "LOOP_HEADER",
    "DATA_BLOCK_HEADING", "DIGIT", "EXPONENT", "CHAR_STRING", "SEMI_COLON_TEXT_FIELD",
    "EOL", "DOUBLE_QUOTE", "SINGLE_QUOTE", "ORDINARY_CHAR", "NON_BLANK_CHAR_",
    "TEXT_LEAD_CHAR", "ANY_PRINT_CHAR", "DATA_", "LOOP_", "GLOBAL_", "STOP_",
    "VERSION", "NON_BLANK_CHAR", "SINGLE_QUOTED_STRING", "DOUBLE_QUOTED_STRING",
    "UNQUOTED_STRING", "'.'", "'?'", "'-'", "'+'", "'('", "')'"
]




class cifParser(Parser):
    grammarFileName = "cif.g"
    antlr_version = version_str_to_tuple("3.1.2")
    antlr_version_str = "3.1.2"
    tokenNames = tokenNames

    def __init__(self, input, state=None):
        if state is None:
            state = RecognizerSharedState()

        Parser.__init__(self, input, state)


        self.dfa12 = self.DFA12(
            self, 12,
            eot = self.DFA12_eot,
            eof = self.DFA12_eof,
            min = self.DFA12_min,
            max = self.DFA12_max,
            accept = self.DFA12_accept,
            special = self.DFA12_special,
            transition = self.DFA12_transition
            )

        self.dfa21 = self.DFA21(
            self, 21,
            eot = self.DFA21_eot,
            eof = self.DFA21_eof,
            min = self.DFA21_min,
            max = self.DFA21_max,
            accept = self.DFA21_accept,
            special = self.DFA21_special,
            transition = self.DFA21_transition
            )

        self.dfa19 = self.DFA19(
            self, 19,
            eot = self.DFA19_eot,
            eof = self.DFA19_eof,
            min = self.DFA19_min,
            max = self.DFA19_max,
            accept = self.DFA19_accept,
            special = self.DFA19_special,
            transition = self.DFA19_transition
            )

        self.dfa22 = self.DFA22(
            self, 22,
            eot = self.DFA22_eot,
            eof = self.DFA22_eof,
            min = self.DFA22_min,
            max = self.DFA22_max,
            accept = self.DFA22_accept,
            special = self.DFA22_special,
            transition = self.DFA22_transition
            )

        self.dfa24 = self.DFA24(
            self, 24,
            eot = self.DFA24_eot,
            eof = self.DFA24_eof,
            min = self.DFA24_min,
            max = self.DFA24_max,
            accept = self.DFA24_accept,
            special = self.DFA24_special,
            transition = self.DFA24_transition
            )













    # $ANTLR start "parse"
    # cif.g:27:1: parse[builder] : cif ;
    def parse(self, builder):

        self.builder = builder
        try:
            try:
                # cif.g:29:2: ( cif )
                # cif.g:29:4: cif
                pass
                self._state.following.append(self.FOLLOW_cif_in_parse44)
                self.cif()

                self._state.following.pop()




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "parse"


    # $ANTLR start "cif"
    # cif.g:34:1: cif : ( COMMENTS )? ( WHITESPACE )* ( data_block ( ( WHITESPACE )* data_block )* ( WHITESPACE )? )? EOF ;
    def cif(self, ):

        try:
            try:
                # cif.g:35:2: ( ( COMMENTS )? ( WHITESPACE )* ( data_block ( ( WHITESPACE )* data_block )* ( WHITESPACE )? )? EOF )
                # cif.g:35:4: ( COMMENTS )? ( WHITESPACE )* ( data_block ( ( WHITESPACE )* data_block )* ( WHITESPACE )? )? EOF
                pass
                # cif.g:35:4: ( COMMENTS )?
                alt1 = 2
                LA1_0 = self.input.LA(1)

                if (LA1_0 == COMMENTS) :
                    alt1 = 1
                if alt1 == 1:
                    # cif.g:35:5: COMMENTS
                    pass
                    self.match(self.input, COMMENTS, self.FOLLOW_COMMENTS_in_cif57)



                # cif.g:35:16: ( WHITESPACE )*
                while True: #loop2
                    alt2 = 2
                    LA2_0 = self.input.LA(1)

                    if (LA2_0 == WHITESPACE) :
                        alt2 = 1


                    if alt2 == 1:
                        # cif.g:35:17: WHITESPACE
                        pass
                        self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_cif62)


                    else:
                        break #loop2


                # cif.g:35:30: ( data_block ( ( WHITESPACE )* data_block )* ( WHITESPACE )? )?
                alt6 = 2
                LA6_0 = self.input.LA(1)

                if (LA6_0 == DATA_BLOCK_HEADING) :
                    alt6 = 1
                if alt6 == 1:
                    # cif.g:35:32: data_block ( ( WHITESPACE )* data_block )* ( WHITESPACE )?
                    pass
                    self._state.following.append(self.FOLLOW_data_block_in_cif68)
                    self.data_block()

                    self._state.following.pop()
                    # cif.g:35:43: ( ( WHITESPACE )* data_block )*
                    while True: #loop4
                        alt4 = 2
                        LA4_0 = self.input.LA(1)

                        if (LA4_0 == WHITESPACE) :
                            LA4_1 = self.input.LA(2)

                            if (LA4_1 == WHITESPACE or LA4_1 == DATA_BLOCK_HEADING) :
                                alt4 = 1


                        elif (LA4_0 == DATA_BLOCK_HEADING) :
                            alt4 = 1


                        if alt4 == 1:
                            # cif.g:35:45: ( WHITESPACE )* data_block
                            pass
                            # cif.g:35:45: ( WHITESPACE )*
                            while True: #loop3
                                alt3 = 2
                                LA3_0 = self.input.LA(1)

                                if (LA3_0 == WHITESPACE) :
                                    alt3 = 1


                                if alt3 == 1:
                                    # cif.g:35:45: WHITESPACE
                                    pass
                                    self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_cif72)


                                else:
                                    break #loop3


                            self._state.following.append(self.FOLLOW_data_block_in_cif75)
                            self.data_block()

                            self._state.following.pop()


                        else:
                            break #loop4


                    # cif.g:35:71: ( WHITESPACE )?
                    alt5 = 2
                    LA5_0 = self.input.LA(1)

                    if (LA5_0 == WHITESPACE) :
                        alt5 = 1
                    if alt5 == 1:
                        # cif.g:35:72: WHITESPACE
                        pass
                        self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_cif81)






                self.match(self.input, EOF, self.FOLLOW_EOF_in_cif88)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "cif"


    # $ANTLR start "loop_body"
    # cif.g:38:1: loop_body : v1= value ( WHITESPACE v2= value )* ;
    def loop_body(self, ):

        v1 = None

        v2 = None


        self.curr_loop_values = flex.std_string()
        try:
            try:
                # cif.g:40:2: (v1= value ( WHITESPACE v2= value )* )
                # cif.g:40:4: v1= value ( WHITESPACE v2= value )*
                pass
                self._state.following.append(self.FOLLOW_value_in_loop_body106)
                v1 = self.value()

                self._state.following.pop()
                #action start
                self.curr_loop_values.append(str(((v1 is not None) and [self.input.toString(v1.start,v1.stop)] or [None])[0]))
                #action end
                # cif.g:42:8: ( WHITESPACE v2= value )*
                while True: #loop7
                    alt7 = 2
                    LA7_0 = self.input.LA(1)

                    if (LA7_0 == WHITESPACE) :
                        LA7_1 = self.input.LA(2)

                        if (LA7_1 == DIGIT or (CHAR_STRING <= LA7_1 <= SEMI_COLON_TEXT_FIELD) or (31 <= LA7_1 <= 34)) :
                            alt7 = 1




                    if alt7 == 1:
                        # cif.g:42:10: WHITESPACE v2= value
                        pass
                        self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_loop_body119)
                        self._state.following.append(self.FOLLOW_value_in_loop_body132)
                        v2 = self.value()

                        self._state.following.pop()
                        #action start
                        self.curr_loop_values.append(str(((v2 is not None) and [self.input.toString(v2.start,v2.stop)] or [None])[0]))
                        #action end


                    else:
                        break #loop7






            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "loop_body"


    # $ANTLR start "save_frame"
    # cif.g:48:1: save_frame : SAVE_FRAME_HEADING ( WHITESPACE data_items )+ WHITESPACE SAVE_ ;
    def save_frame(self, ):

        try:
            try:
                # cif.g:49:2: ( SAVE_FRAME_HEADING ( WHITESPACE data_items )+ WHITESPACE SAVE_ )
                # cif.g:49:4: SAVE_FRAME_HEADING ( WHITESPACE data_items )+ WHITESPACE SAVE_
                pass
                self.match(self.input, SAVE_FRAME_HEADING, self.FOLLOW_SAVE_FRAME_HEADING_in_save_frame156)
                # cif.g:49:23: ( WHITESPACE data_items )+
                cnt8 = 0
                while True: #loop8
                    alt8 = 2
                    LA8_0 = self.input.LA(1)

                    if (LA8_0 == WHITESPACE) :
                        LA8_1 = self.input.LA(2)

                        if ((TAG <= LA8_1 <= LOOP_HEADER)) :
                            alt8 = 1




                    if alt8 == 1:
                        # cif.g:49:25: WHITESPACE data_items
                        pass
                        self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_save_frame160)
                        self._state.following.append(self.FOLLOW_data_items_in_save_frame162)
                        self.data_items()

                        self._state.following.pop()


                    else:
                        if cnt8 >= 1:
                            break #loop8

                        eee = EarlyExitException(8, self.input)
                        raise eee

                    cnt8 += 1


                self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_save_frame167)
                self.match(self.input, SAVE_, self.FOLLOW_SAVE__in_save_frame169)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "save_frame"


    # $ANTLR start "data_items"
    # cif.g:51:1: data_items : ( TAG WHITESPACE value | LOOP_HEADER loop_body );
    def data_items(self, ):

        TAG1 = None
        LOOP_HEADER3 = None
        value2 = None


        try:
            try:
                # cif.g:52:2: ( TAG WHITESPACE value | LOOP_HEADER loop_body )
                alt9 = 2
                LA9_0 = self.input.LA(1)

                if (LA9_0 == TAG) :
                    alt9 = 1
                elif (LA9_0 == LOOP_HEADER) :
                    alt9 = 2
                else:
                    nvae = NoViableAltException("", 9, 0, self.input)

                    raise nvae

                if alt9 == 1:
                    # cif.g:52:4: TAG WHITESPACE value
                    pass
                    TAG1=self.match(self.input, TAG, self.FOLLOW_TAG_in_data_items179)
                    self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_data_items181)
                    self._state.following.append(self.FOLLOW_value_in_data_items183)
                    value2 = self.value()

                    self._state.following.pop()
                    #action start
                    self.builder.add_data_item(TAG1.text, ((value2 is not None) and [self.input.toString(value2.start,value2.stop)] or [None])[0])
                    #action end


                elif alt9 == 2:
                    # cif.g:54:10: LOOP_HEADER loop_body
                    pass
                    LOOP_HEADER3=self.match(self.input, LOOP_HEADER, self.FOLLOW_LOOP_HEADER_in_data_items196)
                    self._state.following.append(self.FOLLOW_loop_body_in_data_items198)
                    self.loop_body()

                    self._state.following.pop()
                    #action start

                    self.builder.add_loop(LOOP_HEADER3.text, data=self.curr_loop_values)

                    #action end



            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "data_items"


    # $ANTLR start "data_block"
    # cif.g:60:1: data_block : DATA_BLOCK_HEADING ( ( WHITESPACE )+ ( data_items | save_frame ) )* ;
    def data_block(self, ):

        DATA_BLOCK_HEADING4 = None

        try:
            try:
                # cif.g:61:2: ( DATA_BLOCK_HEADING ( ( WHITESPACE )+ ( data_items | save_frame ) )* )
                # cif.g:61:4: DATA_BLOCK_HEADING ( ( WHITESPACE )+ ( data_items | save_frame ) )*
                pass
                DATA_BLOCK_HEADING4=self.match(self.input, DATA_BLOCK_HEADING, self.FOLLOW_DATA_BLOCK_HEADING_in_data_block211)
                #action start
                self.builder.add_data_block(DATA_BLOCK_HEADING4.text)
                #action end
                # cif.g:63:8: ( ( WHITESPACE )+ ( data_items | save_frame ) )*
                while True: #loop12
                    alt12 = 2
                    alt12 = self.dfa12.predict(self.input)
                    if alt12 == 1:
                        # cif.g:63:10: ( WHITESPACE )+ ( data_items | save_frame )
                        pass
                        # cif.g:63:10: ( WHITESPACE )+
                        cnt10 = 0
                        while True: #loop10
                            alt10 = 2
                            LA10_0 = self.input.LA(1)

                            if (LA10_0 == WHITESPACE) :
                                alt10 = 1


                            if alt10 == 1:
                                # cif.g:63:10: WHITESPACE
                                pass
                                self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_data_block224)


                            else:
                                if cnt10 >= 1:
                                    break #loop10

                                eee = EarlyExitException(10, self.input)
                                raise eee

                            cnt10 += 1


                        # cif.g:63:22: ( data_items | save_frame )
                        alt11 = 2
                        LA11_0 = self.input.LA(1)

                        if ((TAG <= LA11_0 <= LOOP_HEADER)) :
                            alt11 = 1
                        elif (LA11_0 == SAVE_FRAME_HEADING) :
                            alt11 = 2
                        else:
                            nvae = NoViableAltException("", 11, 0, self.input)

                            raise nvae

                        if alt11 == 1:
                            # cif.g:63:24: data_items
                            pass
                            self._state.following.append(self.FOLLOW_data_items_in_data_block229)
                            self.data_items()

                            self._state.following.pop()


                        elif alt11 == 2:
                            # cif.g:63:37: save_frame
                            pass
                            self._state.following.append(self.FOLLOW_save_frame_in_data_block233)
                            self.save_frame()

                            self._state.following.pop()





                    else:
                        break #loop12






            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "data_block"


    # $ANTLR start "inapplicable"
    # cif.g:70:1: inapplicable : '.' ;
    def inapplicable(self, ):

        try:
            try:
                # cif.g:71:2: ( '.' )
                # cif.g:71:4: '.'
                pass
                self.match(self.input, 31, self.FOLLOW_31_in_inapplicable252)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "inapplicable"


    # $ANTLR start "unknown"
    # cif.g:73:1: unknown : '?' ;
    def unknown(self, ):

        try:
            try:
                # cif.g:73:9: ( '?' )
                # cif.g:73:11: '?'
                pass
                self.match(self.input, 32, self.FOLLOW_32_in_unknown261)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "unknown"

    class value_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)





    # $ANTLR start "value"
    # cif.g:75:1: value : ( inapplicable | unknown | '-' | char_string | numeric | text_field );
    def value(self, ):

        retval = self.value_return()
        retval.start = self.input.LT(1)

        try:
            try:
                # cif.g:75:8: ( inapplicable | unknown | '-' | char_string | numeric | text_field )
                alt13 = 6
                LA13 = self.input.LA(1)
                if LA13 == 31:
                    LA13_1 = self.input.LA(2)

                    if (LA13_1 == DIGIT) :
                        alt13 = 5
                    elif (LA13_1 == EOF or LA13_1 == WHITESPACE or LA13_1 == DATA_BLOCK_HEADING) :
                        alt13 = 1
                    else:
                        nvae = NoViableAltException("", 13, 1, self.input)

                        raise nvae

                elif LA13 == 32:
                    alt13 = 2
                elif LA13 == 33:
                    LA13_3 = self.input.LA(2)

                    if (LA13_3 == EOF or LA13_3 == WHITESPACE or LA13_3 == DATA_BLOCK_HEADING) :
                        alt13 = 3
                    elif (LA13_3 == DIGIT or LA13_3 == 31) :
                        alt13 = 5
                    else:
                        nvae = NoViableAltException("", 13, 3, self.input)

                        raise nvae

                elif LA13 == CHAR_STRING:
                    alt13 = 4
                elif LA13 == DIGIT or LA13 == 34:
                    alt13 = 5
                elif LA13 == SEMI_COLON_TEXT_FIELD:
                    alt13 = 6
                else:
                    nvae = NoViableAltException("", 13, 0, self.input)

                    raise nvae

                if alt13 == 1:
                    # cif.g:75:10: inapplicable
                    pass
                    self._state.following.append(self.FOLLOW_inapplicable_in_value271)
                    self.inapplicable()

                    self._state.following.pop()


                elif alt13 == 2:
                    # cif.g:75:25: unknown
                    pass
                    self._state.following.append(self.FOLLOW_unknown_in_value275)
                    self.unknown()

                    self._state.following.pop()


                elif alt13 == 3:
                    # cif.g:75:35: '-'
                    pass
                    self.match(self.input, 33, self.FOLLOW_33_in_value279)


                elif alt13 == 4:
                    # cif.g:75:41: char_string
                    pass
                    self._state.following.append(self.FOLLOW_char_string_in_value283)
                    self.char_string()

                    self._state.following.pop()


                elif alt13 == 5:
                    # cif.g:75:56: numeric
                    pass
                    self._state.following.append(self.FOLLOW_numeric_in_value288)
                    self.numeric()

                    self._state.following.pop()


                elif alt13 == 6:
                    # cif.g:75:65: text_field
                    pass
                    self._state.following.append(self.FOLLOW_text_field_in_value291)
                    self.text_field()

                    self._state.following.pop()


                retval.stop = self.input.LT(-1)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return retval

    # $ANTLR end "value"


    # $ANTLR start "unsigned_integer"
    # cif.g:77:1: unsigned_integer : ( DIGIT )+ ;
    def unsigned_integer(self, ):

        try:
            try:
                # cif.g:78:2: ( ( DIGIT )+ )
                # cif.g:78:4: ( DIGIT )+
                pass
                # cif.g:78:4: ( DIGIT )+
                cnt14 = 0
                while True: #loop14
                    alt14 = 2
                    LA14_0 = self.input.LA(1)

                    if (LA14_0 == DIGIT) :
                        alt14 = 1


                    if alt14 == 1:
                        # cif.g:78:5: DIGIT
                        pass
                        self.match(self.input, DIGIT, self.FOLLOW_DIGIT_in_unsigned_integer302)


                    else:
                        if cnt14 >= 1:
                            break #loop14

                        eee = EarlyExitException(14, self.input)
                        raise eee

                    cnt14 += 1






            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "unsigned_integer"


    # $ANTLR start "integer"
    # cif.g:80:1: integer : ( '+' | '-' )? unsigned_integer ;
    def integer(self, ):

        try:
            try:
                # cif.g:80:9: ( ( '+' | '-' )? unsigned_integer )
                # cif.g:80:12: ( '+' | '-' )? unsigned_integer
                pass
                # cif.g:80:12: ( '+' | '-' )?
                alt15 = 2
                LA15_0 = self.input.LA(1)

                if ((33 <= LA15_0 <= 34)) :
                    alt15 = 1
                if alt15 == 1:
                    # cif.g:
                    pass
                    if (33 <= self.input.LA(1) <= 34):
                        self.input.consume()
                        self._state.errorRecovery = False

                    else:
                        mse = MismatchedSetException(None, self.input)
                        raise mse





                self._state.following.append(self.FOLLOW_unsigned_integer_in_integer325)
                self.unsigned_integer()

                self._state.following.pop()




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "integer"


    # $ANTLR start "float_"
    # cif.g:82:1: float_ : ( integer EXPONENT | ( ( '+' | '-' )? ( ( DIGIT )* '.' unsigned_integer ) | ( DIGIT )+ '.' ) ( EXPONENT )? );
    def float_(self, ):

        try:
            try:
                # cif.g:82:8: ( integer EXPONENT | ( ( '+' | '-' )? ( ( DIGIT )* '.' unsigned_integer ) | ( DIGIT )+ '.' ) ( EXPONENT )? )
                alt21 = 2
                alt21 = self.dfa21.predict(self.input)
                if alt21 == 1:
                    # cif.g:82:11: integer EXPONENT
                    pass
                    self._state.following.append(self.FOLLOW_integer_in_float_335)
                    self.integer()

                    self._state.following.pop()
                    self.match(self.input, EXPONENT, self.FOLLOW_EXPONENT_in_float_337)


                elif alt21 == 2:
                    # cif.g:82:30: ( ( '+' | '-' )? ( ( DIGIT )* '.' unsigned_integer ) | ( DIGIT )+ '.' ) ( EXPONENT )?
                    pass
                    # cif.g:82:30: ( ( '+' | '-' )? ( ( DIGIT )* '.' unsigned_integer ) | ( DIGIT )+ '.' )
                    alt19 = 2
                    alt19 = self.dfa19.predict(self.input)
                    if alt19 == 1:
                        # cif.g:82:32: ( '+' | '-' )? ( ( DIGIT )* '.' unsigned_integer )
                        pass
                        # cif.g:82:32: ( '+' | '-' )?
                        alt16 = 2
                        LA16_0 = self.input.LA(1)

                        if ((33 <= LA16_0 <= 34)) :
                            alt16 = 1
                        if alt16 == 1:
                            # cif.g:
                            pass
                            if (33 <= self.input.LA(1) <= 34):
                                self.input.consume()
                                self._state.errorRecovery = False

                            else:
                                mse = MismatchedSetException(None, self.input)
                                raise mse





                        # cif.g:82:47: ( ( DIGIT )* '.' unsigned_integer )
                        # cif.g:82:49: ( DIGIT )* '.' unsigned_integer
                        pass
                        # cif.g:82:49: ( DIGIT )*
                        while True: #loop17
                            alt17 = 2
                            LA17_0 = self.input.LA(1)

                            if (LA17_0 == DIGIT) :
                                alt17 = 1


                            if alt17 == 1:
                                # cif.g:82:50: DIGIT
                                pass
                                self.match(self.input, DIGIT, self.FOLLOW_DIGIT_in_float_357)


                            else:
                                break #loop17


                        self.match(self.input, 31, self.FOLLOW_31_in_float_361)
                        self._state.following.append(self.FOLLOW_unsigned_integer_in_float_363)
                        self.unsigned_integer()

                        self._state.following.pop()





                    elif alt19 == 2:
                        # cif.g:82:82: ( DIGIT )+ '.'
                        pass
                        # cif.g:82:82: ( DIGIT )+
                        cnt18 = 0
                        while True: #loop18
                            alt18 = 2
                            LA18_0 = self.input.LA(1)

                            if (LA18_0 == DIGIT) :
                                alt18 = 1


                            if alt18 == 1:
                                # cif.g:82:83: DIGIT
                                pass
                                self.match(self.input, DIGIT, self.FOLLOW_DIGIT_in_float_369)


                            else:
                                if cnt18 >= 1:
                                    break #loop18

                                eee = EarlyExitException(18, self.input)
                                raise eee

                            cnt18 += 1


                        self.match(self.input, 31, self.FOLLOW_31_in_float_373)



                    # cif.g:82:97: ( EXPONENT )?
                    alt20 = 2
                    LA20_0 = self.input.LA(1)

                    if (LA20_0 == EXPONENT) :
                        alt20 = 1
                    if alt20 == 1:
                        # cif.g:82:98: EXPONENT
                        pass
                        self.match(self.input, EXPONENT, self.FOLLOW_EXPONENT_in_float_378)






            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "float_"


    # $ANTLR start "number"
    # cif.g:84:1: number : ( integer | float_ );
    def number(self, ):

        try:
            try:
                # cif.g:84:9: ( integer | float_ )
                alt22 = 2
                alt22 = self.dfa22.predict(self.input)
                if alt22 == 1:
                    # cif.g:84:11: integer
                    pass
                    self._state.following.append(self.FOLLOW_integer_in_number390)
                    self.integer()

                    self._state.following.pop()


                elif alt22 == 2:
                    # cif.g:84:21: float_
                    pass
                    self._state.following.append(self.FOLLOW_float__in_number394)
                    self.float_()

                    self._state.following.pop()



            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "number"


    # $ANTLR start "numeric"
    # cif.g:86:1: numeric : ( number | ( number '(' ( DIGIT )+ ')' ) );
    def numeric(self, ):

        try:
            try:
                # cif.g:86:9: ( number | ( number '(' ( DIGIT )+ ')' ) )
                alt24 = 2
                alt24 = self.dfa24.predict(self.input)
                if alt24 == 1:
                    # cif.g:86:11: number
                    pass
                    self._state.following.append(self.FOLLOW_number_in_numeric403)
                    self.number()

                    self._state.following.pop()


                elif alt24 == 2:
                    # cif.g:86:20: ( number '(' ( DIGIT )+ ')' )
                    pass
                    # cif.g:86:20: ( number '(' ( DIGIT )+ ')' )
                    # cif.g:86:22: number '(' ( DIGIT )+ ')'
                    pass
                    self._state.following.append(self.FOLLOW_number_in_numeric409)
                    self.number()

                    self._state.following.pop()
                    self.match(self.input, 35, self.FOLLOW_35_in_numeric411)
                    # cif.g:86:33: ( DIGIT )+
                    cnt23 = 0
                    while True: #loop23
                        alt23 = 2
                        LA23_0 = self.input.LA(1)

                        if (LA23_0 == DIGIT) :
                            alt23 = 1


                        if alt23 == 1:
                            # cif.g:86:34: DIGIT
                            pass
                            self.match(self.input, DIGIT, self.FOLLOW_DIGIT_in_numeric414)


                        else:
                            if cnt23 >= 1:
                                break #loop23

                            eee = EarlyExitException(23, self.input)
                            raise eee

                        cnt23 += 1


                    self.match(self.input, 36, self.FOLLOW_36_in_numeric418)






            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "numeric"


    # $ANTLR start "char_string"
    # cif.g:88:1: char_string : CHAR_STRING ;
    def char_string(self, ):

        try:
            try:
                # cif.g:89:2: ( CHAR_STRING )
                # cif.g:89:4: CHAR_STRING
                pass
                self.match(self.input, CHAR_STRING, self.FOLLOW_CHAR_STRING_in_char_string430)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "char_string"


    # $ANTLR start "text_field"
    # cif.g:91:1: text_field : SEMI_COLON_TEXT_FIELD ;
    def text_field(self, ):

        try:
            try:
                # cif.g:92:2: ( SEMI_COLON_TEXT_FIELD )
                # cif.g:92:4: SEMI_COLON_TEXT_FIELD
                pass
                self.match(self.input, SEMI_COLON_TEXT_FIELD, self.FOLLOW_SEMI_COLON_TEXT_FIELD_in_text_field440)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "text_field"


    # Delegated rules


    # lookup tables for DFA #12

    DFA12_eot = DFA.unpack(
        u"\5\uffff"
        )

    DFA12_eof = DFA.unpack(
        u"\2\2\3\uffff"
        )

    DFA12_min = DFA.unpack(
        u"\2\5\1\uffff\1\5\1\uffff"
        )

    DFA12_max = DFA.unpack(
        u"\2\12\1\uffff\1\12\1\uffff"
        )

    DFA12_accept = DFA.unpack(
        u"\2\uffff\1\2\1\uffff\1\1"
        )

    DFA12_special = DFA.unpack(
        u"\5\uffff"
        )


    DFA12_transition = [
        DFA.unpack(u"\1\1\4\uffff\1\2"),
        DFA.unpack(u"\1\3\1\4\1\uffff\2\4\1\2"),
        DFA.unpack(u""),
        DFA.unpack(u"\1\3\1\4\1\uffff\2\4\1\2"),
        DFA.unpack(u"")
    ]

    # class definition for DFA #12

    DFA12 = DFA
    # lookup tables for DFA #21

    DFA21_eot = DFA.unpack(
        u"\6\uffff"
        )

    DFA21_eof = DFA.unpack(
        u"\6\uffff"
        )

    DFA21_min = DFA.unpack(
        u"\3\13\1\uffff\1\13\1\uffff"
        )

    DFA21_max = DFA.unpack(
        u"\1\42\2\37\1\uffff\1\37\1\uffff"
        )

    DFA21_accept = DFA.unpack(
        u"\3\uffff\1\2\1\uffff\1\1"
        )

    DFA21_special = DFA.unpack(
        u"\6\uffff"
        )


    DFA21_transition = [
        DFA.unpack(u"\1\2\23\uffff\1\3\1\uffff\2\1"),
        DFA.unpack(u"\1\4\23\uffff\1\3"),
        DFA.unpack(u"\1\2\1\5\22\uffff\1\3"),
        DFA.unpack(u""),
        DFA.unpack(u"\1\4\1\5\22\uffff\1\3"),
        DFA.unpack(u"")
    ]

    # class definition for DFA #21

    DFA21 = DFA
    # lookup tables for DFA #19

    DFA19_eot = DFA.unpack(
        u"\5\uffff"
        )

    DFA19_eof = DFA.unpack(
        u"\3\uffff\1\4\1\uffff"
        )

    DFA19_min = DFA.unpack(
        u"\1\13\1\uffff\1\13\1\5\1\uffff"
        )

    DFA19_max = DFA.unpack(
        u"\1\42\1\uffff\1\37\1\43\1\uffff"
        )

    DFA19_accept = DFA.unpack(
        u"\1\uffff\1\1\2\uffff\1\2"
        )

    DFA19_special = DFA.unpack(
        u"\5\uffff"
        )


    DFA19_transition = [
        DFA.unpack(u"\1\2\23\uffff\1\1\1\uffff\2\1"),
        DFA.unpack(u""),
        DFA.unpack(u"\1\2\23\uffff\1\3"),
        DFA.unpack(u"\1\4\4\uffff\1\4\1\1\1\4\26\uffff\1\4"),
        DFA.unpack(u"")
    ]

    # class definition for DFA #19

    DFA19 = DFA
    # lookup tables for DFA #22

    DFA22_eot = DFA.unpack(
        u"\6\uffff"
        )

    DFA22_eof = DFA.unpack(
        u"\2\uffff\1\5\1\uffff\1\5\1\uffff"
        )

    DFA22_min = DFA.unpack(
        u"\2\13\1\5\1\uffff\1\5\1\uffff"
        )

    DFA22_max = DFA.unpack(
        u"\1\42\1\37\1\43\1\uffff\1\43\1\uffff"
        )

    DFA22_accept = DFA.unpack(
        u"\3\uffff\1\2\1\uffff\1\1"
        )

    DFA22_special = DFA.unpack(
        u"\6\uffff"
        )


    DFA22_transition = [
        DFA.unpack(u"\1\2\23\uffff\1\3\1\uffff\2\1"),
        DFA.unpack(u"\1\4\23\uffff\1\3"),
        DFA.unpack(u"\1\5\4\uffff\1\5\1\2\1\3\22\uffff\1\3\3\uffff\1\5"),
        DFA.unpack(u""),
        DFA.unpack(u"\1\5\4\uffff\1\5\1\4\1\3\22\uffff\1\3\3\uffff\1\5"),
        DFA.unpack(u"")
    ]

    # class definition for DFA #22

    DFA22 = DFA
    # lookup tables for DFA #24

    DFA24_eot = DFA.unpack(
        u"\13\uffff"
        )

    DFA24_eof = DFA.unpack(
        u"\2\uffff\1\10\1\uffff\3\10\2\uffff\2\10"
        )

    DFA24_min = DFA.unpack(
        u"\2\13\1\5\1\13\3\5\2\uffff\2\5"
        )

    DFA24_max = DFA.unpack(
        u"\1\42\1\37\1\43\1\13\3\43\2\uffff\2\43"
        )

    DFA24_accept = DFA.unpack(
        u"\7\uffff\1\2\1\1\2\uffff"
        )

    DFA24_special = DFA.unpack(
        u"\13\uffff"
        )


    DFA24_transition = [
        DFA.unpack(u"\1\2\23\uffff\1\3\1\uffff\2\1"),
        DFA.unpack(u"\1\4\23\uffff\1\3"),
        DFA.unpack(u"\1\10\4\uffff\1\10\1\2\1\6\22\uffff\1\5\3\uffff\1"
        u"\7"),
        DFA.unpack(u"\1\11"),
        DFA.unpack(u"\1\10\4\uffff\1\10\1\4\1\6\22\uffff\1\3\3\uffff\1"
        u"\7"),
        DFA.unpack(u"\1\10\4\uffff\1\10\1\11\1\12\26\uffff\1\7"),
        DFA.unpack(u"\1\10\4\uffff\1\10\30\uffff\1\7"),
        DFA.unpack(u""),
        DFA.unpack(u""),
        DFA.unpack(u"\1\10\4\uffff\1\10\1\11\1\12\26\uffff\1\7"),
        DFA.unpack(u"\1\10\4\uffff\1\10\30\uffff\1\7")
    ]

    # class definition for DFA #24

    DFA24 = DFA


    FOLLOW_cif_in_parse44 = frozenset([1])
    FOLLOW_COMMENTS_in_cif57 = frozenset([5, 10])
    FOLLOW_WHITESPACE_in_cif62 = frozenset([5, 10])
    FOLLOW_data_block_in_cif68 = frozenset([5, 10])
    FOLLOW_WHITESPACE_in_cif72 = frozenset([5, 10])
    FOLLOW_data_block_in_cif75 = frozenset([5, 10])
    FOLLOW_WHITESPACE_in_cif81 = frozenset([])
    FOLLOW_EOF_in_cif88 = frozenset([1])
    FOLLOW_value_in_loop_body106 = frozenset([1, 5])
    FOLLOW_WHITESPACE_in_loop_body119 = frozenset([11, 13, 14, 31, 32, 33, 34])
    FOLLOW_value_in_loop_body132 = frozenset([1, 5])
    FOLLOW_SAVE_FRAME_HEADING_in_save_frame156 = frozenset([5])
    FOLLOW_WHITESPACE_in_save_frame160 = frozenset([8, 9])
    FOLLOW_data_items_in_save_frame162 = frozenset([5])
    FOLLOW_WHITESPACE_in_save_frame167 = frozenset([7])
    FOLLOW_SAVE__in_save_frame169 = frozenset([1])
    FOLLOW_TAG_in_data_items179 = frozenset([5])
    FOLLOW_WHITESPACE_in_data_items181 = frozenset([11, 13, 14, 31, 32, 33, 34])
    FOLLOW_value_in_data_items183 = frozenset([1])
    FOLLOW_LOOP_HEADER_in_data_items196 = frozenset([11, 13, 14, 31, 32, 33, 34])
    FOLLOW_loop_body_in_data_items198 = frozenset([1])
    FOLLOW_DATA_BLOCK_HEADING_in_data_block211 = frozenset([1, 5])
    FOLLOW_WHITESPACE_in_data_block224 = frozenset([5, 6, 8, 9])
    FOLLOW_data_items_in_data_block229 = frozenset([1, 5])
    FOLLOW_save_frame_in_data_block233 = frozenset([1, 5])
    FOLLOW_31_in_inapplicable252 = frozenset([1])
    FOLLOW_32_in_unknown261 = frozenset([1])
    FOLLOW_inapplicable_in_value271 = frozenset([1])
    FOLLOW_unknown_in_value275 = frozenset([1])
    FOLLOW_33_in_value279 = frozenset([1])
    FOLLOW_char_string_in_value283 = frozenset([1])
    FOLLOW_numeric_in_value288 = frozenset([1])
    FOLLOW_text_field_in_value291 = frozenset([1])
    FOLLOW_DIGIT_in_unsigned_integer302 = frozenset([1, 11])
    FOLLOW_set_in_integer314 = frozenset([11, 33, 34])
    FOLLOW_unsigned_integer_in_integer325 = frozenset([1])
    FOLLOW_integer_in_float_335 = frozenset([12])
    FOLLOW_EXPONENT_in_float_337 = frozenset([1])
    FOLLOW_set_in_float_343 = frozenset([11, 31])
    FOLLOW_DIGIT_in_float_357 = frozenset([11, 31])
    FOLLOW_31_in_float_361 = frozenset([11, 33, 34])
    FOLLOW_unsigned_integer_in_float_363 = frozenset([1, 12])
    FOLLOW_DIGIT_in_float_369 = frozenset([11, 31])
    FOLLOW_31_in_float_373 = frozenset([1, 12])
    FOLLOW_EXPONENT_in_float_378 = frozenset([1])
    FOLLOW_integer_in_number390 = frozenset([1])
    FOLLOW_float__in_number394 = frozenset([1])
    FOLLOW_number_in_numeric403 = frozenset([1])
    FOLLOW_number_in_numeric409 = frozenset([35])
    FOLLOW_35_in_numeric411 = frozenset([11])
    FOLLOW_DIGIT_in_numeric414 = frozenset([11, 36])
    FOLLOW_36_in_numeric418 = frozenset([1])
    FOLLOW_CHAR_STRING_in_char_string430 = frozenset([1])
    FOLLOW_SEMI_COLON_TEXT_FIELD_in_text_field440 = frozenset([1])



def main(argv, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr):
    from antlr3.main import ParserMain
    main = ParserMain("cifLexer", cifParser)
    main.stdin = stdin
    main.stdout = stdout
    main.stderr = stderr
    main.execute(argv)


if __name__ == '__main__':
    main(sys.argv)
