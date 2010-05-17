# $ANTLR 3.1.2 cifProto.g 2010-05-17 16:13:59

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
SEMI_COLON_TEXT_FIELD=14
SINGLE_QUOTED_STRING=28
DOUBLE_QUOTE=16
GLOBAL_=24
ORDINARY_CHAR=18
WHITESPACE=5
SAVE=7
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
SAVE_=23
ANY_PRINT_CHAR=21
TEXT_LEAD_CHAR=20
LOOP_=10
DIGIT=11
UNQUOTED_STRING=30
DATA_=22
DATA_BLOCK_HEADING=9

# token names
tokenNames = [
    "<invalid>", "<EOR>", "<DOWN>", "<UP>",
    "COMMENTS", "WHITESPACE", "SAVE_FRAME_HEADING", "SAVE", "TAG", "DATA_BLOCK_HEADING",
    "LOOP_", "DIGIT", "EXPONENT", "CHAR_STRING", "SEMI_COLON_TEXT_FIELD",
    "EOL", "DOUBLE_QUOTE", "SINGLE_QUOTE", "ORDINARY_CHAR", "NON_BLANK_CHAR_",
    "TEXT_LEAD_CHAR", "ANY_PRINT_CHAR", "DATA_", "SAVE_", "GLOBAL_", "STOP_",
    "VERSION", "NON_BLANK_CHAR", "SINGLE_QUOTED_STRING", "DOUBLE_QUOTED_STRING",
    "UNQUOTED_STRING", "'.'", "'?'", "'-'", "'+'", "'('", "')'"
]




class cifProtoParser(Parser):
    grammarFileName = "cifProto.g"
    antlr_version = version_str_to_tuple("3.1.2")
    antlr_version_str = "3.1.2"
    tokenNames = tokenNames

    def __init__(self, input, state=None):
        if state is None:
            state = RecognizerSharedState()

        Parser.__init__(self, input, state)


        self.dfa8 = self.DFA8(
            self, 8,
            eot = self.DFA8_eot,
            eof = self.DFA8_eof,
            min = self.DFA8_min,
            max = self.DFA8_max,
            accept = self.DFA8_accept,
            special = self.DFA8_special,
            transition = self.DFA8_transition
            )

        self.dfa10 = self.DFA10(
            self, 10,
            eot = self.DFA10_eot,
            eof = self.DFA10_eof,
            min = self.DFA10_min,
            max = self.DFA10_max,
            accept = self.DFA10_accept,
            special = self.DFA10_special,
            transition = self.DFA10_transition
            )

        self.dfa16 = self.DFA16(
            self, 16,
            eot = self.DFA16_eot,
            eof = self.DFA16_eof,
            min = self.DFA16_min,
            max = self.DFA16_max,
            accept = self.DFA16_accept,
            special = self.DFA16_special,
            transition = self.DFA16_transition
            )

        self.dfa18 = self.DFA18(
            self, 18,
            eot = self.DFA18_eot,
            eof = self.DFA18_eof,
            min = self.DFA18_min,
            max = self.DFA18_max,
            accept = self.DFA18_accept,
            special = self.DFA18_special,
            transition = self.DFA18_transition
            )

        self.dfa27 = self.DFA27(
            self, 27,
            eot = self.DFA27_eot,
            eof = self.DFA27_eof,
            min = self.DFA27_min,
            max = self.DFA27_max,
            accept = self.DFA27_accept,
            special = self.DFA27_special,
            transition = self.DFA27_transition
            )

        self.dfa25 = self.DFA25(
            self, 25,
            eot = self.DFA25_eot,
            eof = self.DFA25_eof,
            min = self.DFA25_min,
            max = self.DFA25_max,
            accept = self.DFA25_accept,
            special = self.DFA25_special,
            transition = self.DFA25_transition
            )

        self.dfa28 = self.DFA28(
            self, 28,
            eot = self.DFA28_eot,
            eof = self.DFA28_eof,
            min = self.DFA28_min,
            max = self.DFA28_max,
            accept = self.DFA28_accept,
            special = self.DFA28_special,
            transition = self.DFA28_transition
            )

        self.dfa30 = self.DFA30(
            self, 30,
            eot = self.DFA30_eot,
            eof = self.DFA30_eof,
            min = self.DFA30_min,
            max = self.DFA30_max,
            accept = self.DFA30_accept,
            special = self.DFA30_special,
            transition = self.DFA30_transition
            )












    paraphrases = []
    def getErrorMessage(self, e, tokenNames):
      msg = Parser.getErrorMessage(self, e, tokenNames)
      if len(self.paraphrases) > 0:
        paraphrase = self.paraphrases[-1]
        msg += " " + paraphrase
      return msg



    # $ANTLR start "parse"
    # cifProto.g:40:1: parse[builder] : cif EOF ;
    def parse(self, builder):

        self.builder = builder
        try:
            try:
                # cifProto.g:42:2: ( cif EOF )
                # cifProto.g:42:4: cif EOF
                pass
                self._state.following.append(self.FOLLOW_cif_in_parse50)
                self.cif()

                self._state.following.pop()
                self.match(self.input, EOF, self.FOLLOW_EOF_in_parse52)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "parse"


    # $ANTLR start "cif"
    # cifProto.g:47:1: cif : ( COMMENTS )? ( WHITESPACE )* ( data_block ( ( WHITESPACE )* data_block )* ( WHITESPACE )? )? ;
    def cif(self, ):

        try:
            try:
                # cifProto.g:48:2: ( ( COMMENTS )? ( WHITESPACE )* ( data_block ( ( WHITESPACE )* data_block )* ( WHITESPACE )? )? )
                # cifProto.g:48:4: ( COMMENTS )? ( WHITESPACE )* ( data_block ( ( WHITESPACE )* data_block )* ( WHITESPACE )? )?
                pass
                # cifProto.g:48:4: ( COMMENTS )?
                alt1 = 2
                LA1_0 = self.input.LA(1)

                if (LA1_0 == COMMENTS) :
                    alt1 = 1
                if alt1 == 1:
                    # cifProto.g:48:5: COMMENTS
                    pass
                    self.match(self.input, COMMENTS, self.FOLLOW_COMMENTS_in_cif65)



                # cifProto.g:48:16: ( WHITESPACE )*
                while True: #loop2
                    alt2 = 2
                    LA2_0 = self.input.LA(1)

                    if (LA2_0 == WHITESPACE) :
                        alt2 = 1


                    if alt2 == 1:
                        # cifProto.g:48:17: WHITESPACE
                        pass
                        self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_cif70)


                    else:
                        break #loop2


                # cifProto.g:48:30: ( data_block ( ( WHITESPACE )* data_block )* ( WHITESPACE )? )?
                alt6 = 2
                LA6_0 = self.input.LA(1)

                if (LA6_0 == DATA_BLOCK_HEADING) :
                    alt6 = 1
                if alt6 == 1:
                    # cifProto.g:48:32: data_block ( ( WHITESPACE )* data_block )* ( WHITESPACE )?
                    pass
                    self._state.following.append(self.FOLLOW_data_block_in_cif76)
                    self.data_block()

                    self._state.following.pop()
                    # cifProto.g:48:43: ( ( WHITESPACE )* data_block )*
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
                            # cifProto.g:48:45: ( WHITESPACE )* data_block
                            pass
                            # cifProto.g:48:45: ( WHITESPACE )*
                            while True: #loop3
                                alt3 = 2
                                LA3_0 = self.input.LA(1)

                                if (LA3_0 == WHITESPACE) :
                                    alt3 = 1


                                if alt3 == 1:
                                    # cifProto.g:48:45: WHITESPACE
                                    pass
                                    self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_cif80)


                                else:
                                    break #loop3


                            self._state.following.append(self.FOLLOW_data_block_in_cif83)
                            self.data_block()

                            self._state.following.pop()


                        else:
                            break #loop4


                    # cifProto.g:48:71: ( WHITESPACE )?
                    alt5 = 2
                    LA5_0 = self.input.LA(1)

                    if (LA5_0 == WHITESPACE) :
                        alt5 = 1
                    if alt5 == 1:
                        # cifProto.g:48:72: WHITESPACE
                        pass
                        self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_cif89)










            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "cif"


    # $ANTLR start "loop_body"
    # cifProto.g:51:1: loop_body : v1= value ( ( WHITESPACE )+ v2= value )* ;
    def loop_body(self, ):

        v1 = None

        v2 = None


        self.curr_loop_values = flex.std_string()
        try:
            try:
                # cifProto.g:53:2: (v1= value ( ( WHITESPACE )+ v2= value )* )
                # cifProto.g:53:4: v1= value ( ( WHITESPACE )+ v2= value )*
                pass
                self._state.following.append(self.FOLLOW_value_in_loop_body112)
                v1 = self.value()

                self._state.following.pop()
                #action start
                self.curr_loop_values.append(str(((v1 is not None) and [self.input.toString(v1.start,v1.stop)] or [None])[0]))
                #action end
                # cifProto.g:55:8: ( ( WHITESPACE )+ v2= value )*
                while True: #loop8
                    alt8 = 2
                    alt8 = self.dfa8.predict(self.input)
                    if alt8 == 1:
                        # cifProto.g:55:10: ( WHITESPACE )+ v2= value
                        pass
                        # cifProto.g:55:10: ( WHITESPACE )+
                        cnt7 = 0
                        while True: #loop7
                            alt7 = 2
                            LA7_0 = self.input.LA(1)

                            if (LA7_0 == WHITESPACE) :
                                alt7 = 1


                            if alt7 == 1:
                                # cifProto.g:55:10: WHITESPACE
                                pass
                                self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_loop_body125)


                            else:
                                if cnt7 >= 1:
                                    break #loop7

                                eee = EarlyExitException(7, self.input)
                                raise eee

                            cnt7 += 1


                        self._state.following.append(self.FOLLOW_value_in_loop_body139)
                        v2 = self.value()

                        self._state.following.pop()
                        #action start
                        self.curr_loop_values.append(str(((v2 is not None) and [self.input.toString(v2.start,v2.stop)] or [None])[0]))
                        #action end


                    else:
                        break #loop8


                #action start
                self.paraphrases.pop()
                #action end




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "loop_body"


    # $ANTLR start "save_frame"
    # cifProto.g:62:1: save_frame : SAVE_FRAME_HEADING ( ( WHITESPACE )+ data_items )+ ( WHITESPACE )+ SAVE ;
    def save_frame(self, ):

        try:
            try:
                # cifProto.g:63:2: ( SAVE_FRAME_HEADING ( ( WHITESPACE )+ data_items )+ ( WHITESPACE )+ SAVE )
                # cifProto.g:63:4: SAVE_FRAME_HEADING ( ( WHITESPACE )+ data_items )+ ( WHITESPACE )+ SAVE
                pass
                self.match(self.input, SAVE_FRAME_HEADING, self.FOLLOW_SAVE_FRAME_HEADING_in_save_frame174)
                # cifProto.g:63:23: ( ( WHITESPACE )+ data_items )+
                cnt10 = 0
                while True: #loop10
                    alt10 = 2
                    alt10 = self.dfa10.predict(self.input)
                    if alt10 == 1:
                        # cifProto.g:63:25: ( WHITESPACE )+ data_items
                        pass
                        # cifProto.g:63:25: ( WHITESPACE )+
                        cnt9 = 0
                        while True: #loop9
                            alt9 = 2
                            LA9_0 = self.input.LA(1)

                            if (LA9_0 == WHITESPACE) :
                                alt9 = 1


                            if alt9 == 1:
                                # cifProto.g:63:25: WHITESPACE
                                pass
                                self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_save_frame178)


                            else:
                                if cnt9 >= 1:
                                    break #loop9

                                eee = EarlyExitException(9, self.input)
                                raise eee

                            cnt9 += 1


                        self._state.following.append(self.FOLLOW_data_items_in_save_frame181)
                        self.data_items()

                        self._state.following.pop()


                    else:
                        if cnt10 >= 1:
                            break #loop10

                        eee = EarlyExitException(10, self.input)
                        raise eee

                    cnt10 += 1


                # cifProto.g:63:51: ( WHITESPACE )+
                cnt11 = 0
                while True: #loop11
                    alt11 = 2
                    LA11_0 = self.input.LA(1)

                    if (LA11_0 == WHITESPACE) :
                        alt11 = 1


                    if alt11 == 1:
                        # cifProto.g:63:51: WHITESPACE
                        pass
                        self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_save_frame186)


                    else:
                        if cnt11 >= 1:
                            break #loop11

                        eee = EarlyExitException(11, self.input)
                        raise eee

                    cnt11 += 1


                self.match(self.input, SAVE, self.FOLLOW_SAVE_in_save_frame189)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "save_frame"


    # $ANTLR start "data_items"
    # cifProto.g:65:1: data_items : ( ( TAG WHITESPACE value ) | ( loop_header ( WHITESPACE )* loop_body ) );
    def data_items(self, ):

        TAG1 = None
        value2 = None

        loop_header3 = None


        self.paraphrases.append("in data items")
        try:
            try:
                # cifProto.g:68:2: ( ( TAG WHITESPACE value ) | ( loop_header ( WHITESPACE )* loop_body ) )
                alt13 = 2
                LA13_0 = self.input.LA(1)

                if (LA13_0 == TAG) :
                    alt13 = 1
                elif (LA13_0 == LOOP_) :
                    alt13 = 2
                else:
                    nvae = NoViableAltException("", 13, 0, self.input)

                    raise nvae

                if alt13 == 1:
                    # cifProto.g:68:5: ( TAG WHITESPACE value )
                    pass
                    # cifProto.g:68:5: ( TAG WHITESPACE value )
                    # cifProto.g:68:7: TAG WHITESPACE value
                    pass
                    TAG1=self.match(self.input, TAG, self.FOLLOW_TAG_in_data_items213)
                    self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_data_items215)
                    self._state.following.append(self.FOLLOW_value_in_data_items217)
                    value2 = self.value()

                    self._state.following.pop()
                    #action start
                    self.builder.add_data_item(TAG1.text, ((value2 is not None) and [self.input.toString(value2.start,value2.stop)] or [None])[0])
                    #action end





                elif alt13 == 2:
                    # cifProto.g:71:10: ( loop_header ( WHITESPACE )* loop_body )
                    pass
                    # cifProto.g:71:10: ( loop_header ( WHITESPACE )* loop_body )
                    # cifProto.g:71:12: loop_header ( WHITESPACE )* loop_body
                    pass
                    self._state.following.append(self.FOLLOW_loop_header_in_data_items236)
                    loop_header3 = self.loop_header()

                    self._state.following.pop()
                    # cifProto.g:71:24: ( WHITESPACE )*
                    while True: #loop12
                        alt12 = 2
                        LA12_0 = self.input.LA(1)

                        if (LA12_0 == WHITESPACE) :
                            alt12 = 1


                        if alt12 == 1:
                            # cifProto.g:71:24: WHITESPACE
                            pass
                            self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_data_items238)


                        else:
                            break #loop12


                    self._state.following.append(self.FOLLOW_loop_body_in_data_items241)
                    self.loop_body()

                    self._state.following.pop()
                    #action start

                    try:
                      self.builder.add_loop(((loop_header3 is not None) and [self.input.toString(loop_header3.start,loop_header3.stop)] or [None])[0], data=self.curr_loop_values)
                    except AssertionError, e:
                      print e

                    #action end





                #action start
                self.paraphrases.pop()
                #action end

            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "data_items"


    # $ANTLR start "data_block"
    # cifProto.g:81:1: data_block : DATA_BLOCK_HEADING ( ( WHITESPACE )+ ( data_items | save_frame ) )* ;
    def data_block(self, ):

        DATA_BLOCK_HEADING4 = None

        try:
            try:
                # cifProto.g:82:2: ( DATA_BLOCK_HEADING ( ( WHITESPACE )+ ( data_items | save_frame ) )* )
                # cifProto.g:82:4: DATA_BLOCK_HEADING ( ( WHITESPACE )+ ( data_items | save_frame ) )*
                pass
                DATA_BLOCK_HEADING4=self.match(self.input, DATA_BLOCK_HEADING, self.FOLLOW_DATA_BLOCK_HEADING_in_data_block258)
                #action start
                self.builder.add_data_block(DATA_BLOCK_HEADING4.text)
                #action end
                # cifProto.g:84:8: ( ( WHITESPACE )+ ( data_items | save_frame ) )*
                while True: #loop16
                    alt16 = 2
                    alt16 = self.dfa16.predict(self.input)
                    if alt16 == 1:
                        # cifProto.g:84:10: ( WHITESPACE )+ ( data_items | save_frame )
                        pass
                        # cifProto.g:84:10: ( WHITESPACE )+
                        cnt14 = 0
                        while True: #loop14
                            alt14 = 2
                            LA14_0 = self.input.LA(1)

                            if (LA14_0 == WHITESPACE) :
                                alt14 = 1


                            if alt14 == 1:
                                # cifProto.g:84:10: WHITESPACE
                                pass
                                self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_data_block271)


                            else:
                                if cnt14 >= 1:
                                    break #loop14

                                eee = EarlyExitException(14, self.input)
                                raise eee

                            cnt14 += 1


                        # cifProto.g:84:22: ( data_items | save_frame )
                        alt15 = 2
                        LA15_0 = self.input.LA(1)

                        if (LA15_0 == TAG or LA15_0 == LOOP_) :
                            alt15 = 1
                        elif (LA15_0 == SAVE_FRAME_HEADING) :
                            alt15 = 2
                        else:
                            nvae = NoViableAltException("", 15, 0, self.input)

                            raise nvae

                        if alt15 == 1:
                            # cifProto.g:84:24: data_items
                            pass
                            self._state.following.append(self.FOLLOW_data_items_in_data_block276)
                            self.data_items()

                            self._state.following.pop()


                        elif alt15 == 2:
                            # cifProto.g:84:37: save_frame
                            pass
                            self._state.following.append(self.FOLLOW_save_frame_in_data_block280)
                            self.save_frame()

                            self._state.following.pop()





                    else:
                        break #loop16






            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "data_block"

    class loop_header_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)





    # $ANTLR start "loop_header"
    # cifProto.g:87:1: loop_header : LOOP_ ( ( WHITESPACE )+ TAG )+ WHITESPACE ;
    def loop_header(self, ):

        retval = self.loop_header_return()
        retval.start = self.input.LT(1)

        try:
            try:
                # cifProto.g:88:2: ( LOOP_ ( ( WHITESPACE )+ TAG )+ WHITESPACE )
                # cifProto.g:88:4: LOOP_ ( ( WHITESPACE )+ TAG )+ WHITESPACE
                pass
                self.match(self.input, LOOP_, self.FOLLOW_LOOP__in_loop_header296)
                #action start
                self.paraphrases.append("in loop_")
                #action end
                # cifProto.g:89:3: ( ( WHITESPACE )+ TAG )+
                cnt18 = 0
                while True: #loop18
                    alt18 = 2
                    alt18 = self.dfa18.predict(self.input)
                    if alt18 == 1:
                        # cifProto.g:89:5: ( WHITESPACE )+ TAG
                        pass
                        # cifProto.g:89:5: ( WHITESPACE )+
                        cnt17 = 0
                        while True: #loop17
                            alt17 = 2
                            LA17_0 = self.input.LA(1)

                            if (LA17_0 == WHITESPACE) :
                                alt17 = 1


                            if alt17 == 1:
                                # cifProto.g:89:5: WHITESPACE
                                pass
                                self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_loop_header304)


                            else:
                                if cnt17 >= 1:
                                    break #loop17

                                eee = EarlyExitException(17, self.input)
                                raise eee

                            cnt17 += 1


                        self.match(self.input, TAG, self.FOLLOW_TAG_in_loop_header307)


                    else:
                        if cnt18 >= 1:
                            break #loop18

                        eee = EarlyExitException(18, self.input)
                        raise eee

                    cnt18 += 1


                self.match(self.input, WHITESPACE, self.FOLLOW_WHITESPACE_in_loop_header312)



                retval.stop = self.input.LT(-1)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return retval

    # $ANTLR end "loop_header"


    # $ANTLR start "inapplicable"
    # cifProto.g:96:1: inapplicable : '.' ;
    def inapplicable(self, ):

        try:
            try:
                # cifProto.g:97:2: ( '.' )
                # cifProto.g:97:4: '.'
                pass
                self.match(self.input, 31, self.FOLLOW_31_in_inapplicable327)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "inapplicable"


    # $ANTLR start "unknown"
    # cifProto.g:99:1: unknown : '?' ;
    def unknown(self, ):

        try:
            try:
                # cifProto.g:99:9: ( '?' )
                # cifProto.g:99:11: '?'
                pass
                self.match(self.input, 32, self.FOLLOW_32_in_unknown336)




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
    # cifProto.g:101:1: value : ( inapplicable | unknown | '-' | char_string | numeric | text_field );
    def value(self, ):

        retval = self.value_return()
        retval.start = self.input.LT(1)

        self.paraphrases.append("in value")
        try:
            try:
                # cifProto.g:104:3: ( inapplicable | unknown | '-' | char_string | numeric | text_field )
                alt19 = 6
                LA19 = self.input.LA(1)
                if LA19 == 31:
                    LA19_1 = self.input.LA(2)

                    if (LA19_1 == DIGIT) :
                        alt19 = 5
                    elif (LA19_1 == EOF or LA19_1 == WHITESPACE or LA19_1 == DATA_BLOCK_HEADING) :
                        alt19 = 1
                    else:
                        nvae = NoViableAltException("", 19, 1, self.input)

                        raise nvae

                elif LA19 == 32:
                    alt19 = 2
                elif LA19 == 33:
                    LA19_3 = self.input.LA(2)

                    if (LA19_3 == DIGIT or LA19_3 == 31) :
                        alt19 = 5
                    elif (LA19_3 == EOF or LA19_3 == WHITESPACE or LA19_3 == DATA_BLOCK_HEADING) :
                        alt19 = 3
                    else:
                        nvae = NoViableAltException("", 19, 3, self.input)

                        raise nvae

                elif LA19 == CHAR_STRING:
                    alt19 = 4
                elif LA19 == DIGIT or LA19 == 34:
                    alt19 = 5
                elif LA19 == SEMI_COLON_TEXT_FIELD:
                    alt19 = 6
                else:
                    nvae = NoViableAltException("", 19, 0, self.input)

                    raise nvae

                if alt19 == 1:
                    # cifProto.g:104:5: inapplicable
                    pass
                    self._state.following.append(self.FOLLOW_inapplicable_in_value358)
                    self.inapplicable()

                    self._state.following.pop()


                elif alt19 == 2:
                    # cifProto.g:104:20: unknown
                    pass
                    self._state.following.append(self.FOLLOW_unknown_in_value362)
                    self.unknown()

                    self._state.following.pop()


                elif alt19 == 3:
                    # cifProto.g:104:30: '-'
                    pass
                    self.match(self.input, 33, self.FOLLOW_33_in_value366)


                elif alt19 == 4:
                    # cifProto.g:104:36: char_string
                    pass
                    self._state.following.append(self.FOLLOW_char_string_in_value370)
                    self.char_string()

                    self._state.following.pop()


                elif alt19 == 5:
                    # cifProto.g:104:51: numeric
                    pass
                    self._state.following.append(self.FOLLOW_numeric_in_value375)
                    self.numeric()

                    self._state.following.pop()


                elif alt19 == 6:
                    # cifProto.g:104:60: text_field
                    pass
                    self._state.following.append(self.FOLLOW_text_field_in_value378)
                    self.text_field()

                    self._state.following.pop()


                retval.stop = self.input.LT(-1)

                #action start
                self.paraphrases.pop()
                #action end

            except RecognitionException, re:

                raise re


        finally:

            pass

        return retval

    # $ANTLR end "value"


    # $ANTLR start "unsigned_integer"
    # cifProto.g:110:1: unsigned_integer : ( DIGIT )+ ;
    def unsigned_integer(self, ):

        try:
            try:
                # cifProto.g:111:2: ( ( DIGIT )+ )
                # cifProto.g:111:4: ( DIGIT )+
                pass
                # cifProto.g:111:4: ( DIGIT )+
                cnt20 = 0
                while True: #loop20
                    alt20 = 2
                    LA20_0 = self.input.LA(1)

                    if (LA20_0 == DIGIT) :
                        alt20 = 1


                    if alt20 == 1:
                        # cifProto.g:111:5: DIGIT
                        pass
                        self.match(self.input, DIGIT, self.FOLLOW_DIGIT_in_unsigned_integer399)


                    else:
                        if cnt20 >= 1:
                            break #loop20

                        eee = EarlyExitException(20, self.input)
                        raise eee

                    cnt20 += 1






            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "unsigned_integer"


    # $ANTLR start "integer"
    # cifProto.g:113:1: integer : ( '+' | '-' )? unsigned_integer ;
    def integer(self, ):

        try:
            try:
                # cifProto.g:113:9: ( ( '+' | '-' )? unsigned_integer )
                # cifProto.g:113:12: ( '+' | '-' )? unsigned_integer
                pass
                # cifProto.g:113:12: ( '+' | '-' )?
                alt21 = 2
                LA21_0 = self.input.LA(1)

                if ((33 <= LA21_0 <= 34)) :
                    alt21 = 1
                if alt21 == 1:
                    # cifProto.g:
                    pass
                    if (33 <= self.input.LA(1) <= 34):
                        self.input.consume()
                        self._state.errorRecovery = False

                    else:
                        mse = MismatchedSetException(None, self.input)
                        raise mse





                self._state.following.append(self.FOLLOW_unsigned_integer_in_integer422)
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
    # cifProto.g:115:1: float_ : ( integer EXPONENT | ( ( '+' | '-' )? ( ( DIGIT )* '.' unsigned_integer ) | ( DIGIT )+ '.' ) ( EXPONENT )? );
    def float_(self, ):

        try:
            try:
                # cifProto.g:115:8: ( integer EXPONENT | ( ( '+' | '-' )? ( ( DIGIT )* '.' unsigned_integer ) | ( DIGIT )+ '.' ) ( EXPONENT )? )
                alt27 = 2
                alt27 = self.dfa27.predict(self.input)
                if alt27 == 1:
                    # cifProto.g:115:11: integer EXPONENT
                    pass
                    self._state.following.append(self.FOLLOW_integer_in_float_432)
                    self.integer()

                    self._state.following.pop()
                    self.match(self.input, EXPONENT, self.FOLLOW_EXPONENT_in_float_434)


                elif alt27 == 2:
                    # cifProto.g:115:30: ( ( '+' | '-' )? ( ( DIGIT )* '.' unsigned_integer ) | ( DIGIT )+ '.' ) ( EXPONENT )?
                    pass
                    # cifProto.g:115:30: ( ( '+' | '-' )? ( ( DIGIT )* '.' unsigned_integer ) | ( DIGIT )+ '.' )
                    alt25 = 2
                    alt25 = self.dfa25.predict(self.input)
                    if alt25 == 1:
                        # cifProto.g:115:32: ( '+' | '-' )? ( ( DIGIT )* '.' unsigned_integer )
                        pass
                        # cifProto.g:115:32: ( '+' | '-' )?
                        alt22 = 2
                        LA22_0 = self.input.LA(1)

                        if ((33 <= LA22_0 <= 34)) :
                            alt22 = 1
                        if alt22 == 1:
                            # cifProto.g:
                            pass
                            if (33 <= self.input.LA(1) <= 34):
                                self.input.consume()
                                self._state.errorRecovery = False

                            else:
                                mse = MismatchedSetException(None, self.input)
                                raise mse





                        # cifProto.g:115:47: ( ( DIGIT )* '.' unsigned_integer )
                        # cifProto.g:115:49: ( DIGIT )* '.' unsigned_integer
                        pass
                        # cifProto.g:115:49: ( DIGIT )*
                        while True: #loop23
                            alt23 = 2
                            LA23_0 = self.input.LA(1)

                            if (LA23_0 == DIGIT) :
                                alt23 = 1


                            if alt23 == 1:
                                # cifProto.g:115:50: DIGIT
                                pass
                                self.match(self.input, DIGIT, self.FOLLOW_DIGIT_in_float_454)


                            else:
                                break #loop23


                        self.match(self.input, 31, self.FOLLOW_31_in_float_458)
                        self._state.following.append(self.FOLLOW_unsigned_integer_in_float_460)
                        self.unsigned_integer()

                        self._state.following.pop()





                    elif alt25 == 2:
                        # cifProto.g:115:82: ( DIGIT )+ '.'
                        pass
                        # cifProto.g:115:82: ( DIGIT )+
                        cnt24 = 0
                        while True: #loop24
                            alt24 = 2
                            LA24_0 = self.input.LA(1)

                            if (LA24_0 == DIGIT) :
                                alt24 = 1


                            if alt24 == 1:
                                # cifProto.g:115:83: DIGIT
                                pass
                                self.match(self.input, DIGIT, self.FOLLOW_DIGIT_in_float_466)


                            else:
                                if cnt24 >= 1:
                                    break #loop24

                                eee = EarlyExitException(24, self.input)
                                raise eee

                            cnt24 += 1


                        self.match(self.input, 31, self.FOLLOW_31_in_float_470)



                    # cifProto.g:115:97: ( EXPONENT )?
                    alt26 = 2
                    LA26_0 = self.input.LA(1)

                    if (LA26_0 == EXPONENT) :
                        alt26 = 1
                    if alt26 == 1:
                        # cifProto.g:115:98: EXPONENT
                        pass
                        self.match(self.input, EXPONENT, self.FOLLOW_EXPONENT_in_float_475)






            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "float_"


    # $ANTLR start "number"
    # cifProto.g:117:1: number : ( integer | float_ );
    def number(self, ):

        try:
            try:
                # cifProto.g:117:9: ( integer | float_ )
                alt28 = 2
                alt28 = self.dfa28.predict(self.input)
                if alt28 == 1:
                    # cifProto.g:117:11: integer
                    pass
                    self._state.following.append(self.FOLLOW_integer_in_number487)
                    self.integer()

                    self._state.following.pop()


                elif alt28 == 2:
                    # cifProto.g:117:21: float_
                    pass
                    self._state.following.append(self.FOLLOW_float__in_number491)
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
    # cifProto.g:119:1: numeric : ( number | ( number '(' ( DIGIT )+ ')' ) );
    def numeric(self, ):

        try:
            try:
                # cifProto.g:119:9: ( number | ( number '(' ( DIGIT )+ ')' ) )
                alt30 = 2
                alt30 = self.dfa30.predict(self.input)
                if alt30 == 1:
                    # cifProto.g:119:11: number
                    pass
                    self._state.following.append(self.FOLLOW_number_in_numeric500)
                    self.number()

                    self._state.following.pop()


                elif alt30 == 2:
                    # cifProto.g:119:20: ( number '(' ( DIGIT )+ ')' )
                    pass
                    # cifProto.g:119:20: ( number '(' ( DIGIT )+ ')' )
                    # cifProto.g:119:22: number '(' ( DIGIT )+ ')'
                    pass
                    self._state.following.append(self.FOLLOW_number_in_numeric506)
                    self.number()

                    self._state.following.pop()
                    self.match(self.input, 35, self.FOLLOW_35_in_numeric508)
                    # cifProto.g:119:33: ( DIGIT )+
                    cnt29 = 0
                    while True: #loop29
                        alt29 = 2
                        LA29_0 = self.input.LA(1)

                        if (LA29_0 == DIGIT) :
                            alt29 = 1


                        if alt29 == 1:
                            # cifProto.g:119:34: DIGIT
                            pass
                            self.match(self.input, DIGIT, self.FOLLOW_DIGIT_in_numeric511)


                        else:
                            if cnt29 >= 1:
                                break #loop29

                            eee = EarlyExitException(29, self.input)
                            raise eee

                        cnt29 += 1


                    self.match(self.input, 36, self.FOLLOW_36_in_numeric515)






            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "numeric"


    # $ANTLR start "char_string"
    # cifProto.g:121:1: char_string : CHAR_STRING ;
    def char_string(self, ):

        try:
            try:
                # cifProto.g:122:2: ( CHAR_STRING )
                # cifProto.g:122:4: CHAR_STRING
                pass
                self.match(self.input, CHAR_STRING, self.FOLLOW_CHAR_STRING_in_char_string527)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "char_string"


    # $ANTLR start "text_field"
    # cifProto.g:124:1: text_field : SEMI_COLON_TEXT_FIELD ;
    def text_field(self, ):

        try:
            try:
                # cifProto.g:125:2: ( SEMI_COLON_TEXT_FIELD )
                # cifProto.g:125:4: SEMI_COLON_TEXT_FIELD
                pass
                self.match(self.input, SEMI_COLON_TEXT_FIELD, self.FOLLOW_SEMI_COLON_TEXT_FIELD_in_text_field537)




            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
        finally:

            pass

        return

    # $ANTLR end "text_field"


    # Delegated rules


    # lookup tables for DFA #8

    DFA8_eot = DFA.unpack(
        u"\5\uffff"
        )

    DFA8_eof = DFA.unpack(
        u"\2\2\3\uffff"
        )

    DFA8_min = DFA.unpack(
        u"\2\5\1\uffff\1\5\1\uffff"
        )

    DFA8_max = DFA.unpack(
        u"\1\11\1\42\1\uffff\1\42\1\uffff"
        )

    DFA8_accept = DFA.unpack(
        u"\2\uffff\1\2\1\uffff\1\1"
        )

    DFA8_special = DFA.unpack(
        u"\5\uffff"
        )


    DFA8_transition = [
        DFA.unpack(u"\1\1\3\uffff\1\2"),
        DFA.unpack(u"\1\3\5\2\1\4\1\uffff\2\4\20\uffff\4\4"),
        DFA.unpack(u""),
        DFA.unpack(u"\1\3\5\2\1\4\1\uffff\2\4\20\uffff\4\4"),
        DFA.unpack(u"")
    ]

    # class definition for DFA #8

    DFA8 = DFA
    # lookup tables for DFA #10

    DFA10_eot = DFA.unpack(
        u"\4\uffff"
        )

    DFA10_eof = DFA.unpack(
        u"\4\uffff"
        )

    DFA10_min = DFA.unpack(
        u"\2\5\2\uffff"
        )

    DFA10_max = DFA.unpack(
        u"\1\5\1\12\2\uffff"
        )

    DFA10_accept = DFA.unpack(
        u"\2\uffff\1\1\1\2"
        )

    DFA10_special = DFA.unpack(
        u"\4\uffff"
        )


    DFA10_transition = [
        DFA.unpack(u"\1\1"),
        DFA.unpack(u"\1\1\1\uffff\1\3\1\2\1\uffff\1\2"),
        DFA.unpack(u""),
        DFA.unpack(u"")
    ]

    # class definition for DFA #10

    DFA10 = DFA
    # lookup tables for DFA #16

    DFA16_eot = DFA.unpack(
        u"\5\uffff"
        )

    DFA16_eof = DFA.unpack(
        u"\2\2\3\uffff"
        )

    DFA16_min = DFA.unpack(
        u"\2\5\1\uffff\1\5\1\uffff"
        )

    DFA16_max = DFA.unpack(
        u"\1\11\1\12\1\uffff\1\12\1\uffff"
        )

    DFA16_accept = DFA.unpack(
        u"\2\uffff\1\2\1\uffff\1\1"
        )

    DFA16_special = DFA.unpack(
        u"\5\uffff"
        )


    DFA16_transition = [
        DFA.unpack(u"\1\1\3\uffff\1\2"),
        DFA.unpack(u"\1\3\1\4\1\uffff\1\4\1\2\1\4"),
        DFA.unpack(u""),
        DFA.unpack(u"\1\3\1\4\1\uffff\1\4\1\2\1\4"),
        DFA.unpack(u"")
    ]

    # class definition for DFA #16

    DFA16 = DFA
    # lookup tables for DFA #18

    DFA18_eot = DFA.unpack(
        u"\5\uffff"
        )

    DFA18_eof = DFA.unpack(
        u"\5\uffff"
        )

    DFA18_min = DFA.unpack(
        u"\2\5\1\uffff\1\5\1\uffff"
        )

    DFA18_max = DFA.unpack(
        u"\1\5\1\42\1\uffff\1\42\1\uffff"
        )

    DFA18_accept = DFA.unpack(
        u"\2\uffff\1\1\1\uffff\1\2"
        )

    DFA18_special = DFA.unpack(
        u"\5\uffff"
        )


    DFA18_transition = [
        DFA.unpack(u"\1\1"),
        DFA.unpack(u"\1\3\2\uffff\1\2\2\uffff\1\4\1\uffff\2\4\20\uffff"
        u"\4\4"),
        DFA.unpack(u""),
        DFA.unpack(u"\1\3\2\uffff\1\2\2\uffff\1\4\1\uffff\2\4\20\uffff"
        u"\4\4"),
        DFA.unpack(u"")
    ]

    # class definition for DFA #18

    DFA18 = DFA
    # lookup tables for DFA #27

    DFA27_eot = DFA.unpack(
        u"\6\uffff"
        )

    DFA27_eof = DFA.unpack(
        u"\6\uffff"
        )

    DFA27_min = DFA.unpack(
        u"\3\13\1\uffff\1\13\1\uffff"
        )

    DFA27_max = DFA.unpack(
        u"\1\42\2\37\1\uffff\1\37\1\uffff"
        )

    DFA27_accept = DFA.unpack(
        u"\3\uffff\1\2\1\uffff\1\1"
        )

    DFA27_special = DFA.unpack(
        u"\6\uffff"
        )


    DFA27_transition = [
        DFA.unpack(u"\1\2\23\uffff\1\3\1\uffff\2\1"),
        DFA.unpack(u"\1\4\23\uffff\1\3"),
        DFA.unpack(u"\1\2\1\5\22\uffff\1\3"),
        DFA.unpack(u""),
        DFA.unpack(u"\1\4\1\5\22\uffff\1\3"),
        DFA.unpack(u"")
    ]

    # class definition for DFA #27

    DFA27 = DFA
    # lookup tables for DFA #25

    DFA25_eot = DFA.unpack(
        u"\5\uffff"
        )

    DFA25_eof = DFA.unpack(
        u"\3\uffff\1\4\1\uffff"
        )

    DFA25_min = DFA.unpack(
        u"\1\13\1\uffff\1\13\1\5\1\uffff"
        )

    DFA25_max = DFA.unpack(
        u"\1\42\1\uffff\1\37\1\43\1\uffff"
        )

    DFA25_accept = DFA.unpack(
        u"\1\uffff\1\1\2\uffff\1\2"
        )

    DFA25_special = DFA.unpack(
        u"\5\uffff"
        )


    DFA25_transition = [
        DFA.unpack(u"\1\2\23\uffff\1\1\1\uffff\2\1"),
        DFA.unpack(u""),
        DFA.unpack(u"\1\2\23\uffff\1\3"),
        DFA.unpack(u"\1\4\3\uffff\1\4\1\uffff\1\1\1\4\26\uffff\1\4"),
        DFA.unpack(u"")
    ]

    # class definition for DFA #25

    DFA25 = DFA
    # lookup tables for DFA #28

    DFA28_eot = DFA.unpack(
        u"\6\uffff"
        )

    DFA28_eof = DFA.unpack(
        u"\2\uffff\1\5\1\uffff\1\5\1\uffff"
        )

    DFA28_min = DFA.unpack(
        u"\2\13\1\5\1\uffff\1\5\1\uffff"
        )

    DFA28_max = DFA.unpack(
        u"\1\42\1\37\1\43\1\uffff\1\43\1\uffff"
        )

    DFA28_accept = DFA.unpack(
        u"\3\uffff\1\2\1\uffff\1\1"
        )

    DFA28_special = DFA.unpack(
        u"\6\uffff"
        )


    DFA28_transition = [
        DFA.unpack(u"\1\2\23\uffff\1\3\1\uffff\2\1"),
        DFA.unpack(u"\1\4\23\uffff\1\3"),
        DFA.unpack(u"\1\5\3\uffff\1\5\1\uffff\1\2\1\3\22\uffff\1\3\3\uffff"
        u"\1\5"),
        DFA.unpack(u""),
        DFA.unpack(u"\1\5\3\uffff\1\5\1\uffff\1\4\1\3\22\uffff\1\3\3\uffff"
        u"\1\5"),
        DFA.unpack(u"")
    ]

    # class definition for DFA #28

    DFA28 = DFA
    # lookup tables for DFA #30

    DFA30_eot = DFA.unpack(
        u"\13\uffff"
        )

    DFA30_eof = DFA.unpack(
        u"\2\uffff\1\10\1\uffff\3\10\2\uffff\2\10"
        )

    DFA30_min = DFA.unpack(
        u"\2\13\1\5\1\13\3\5\2\uffff\2\5"
        )

    DFA30_max = DFA.unpack(
        u"\1\42\1\37\1\43\1\13\3\43\2\uffff\2\43"
        )

    DFA30_accept = DFA.unpack(
        u"\7\uffff\1\2\1\1\2\uffff"
        )

    DFA30_special = DFA.unpack(
        u"\13\uffff"
        )


    DFA30_transition = [
        DFA.unpack(u"\1\2\23\uffff\1\3\1\uffff\2\1"),
        DFA.unpack(u"\1\4\23\uffff\1\3"),
        DFA.unpack(u"\1\10\3\uffff\1\10\1\uffff\1\2\1\5\22\uffff\1\6\3"
        u"\uffff\1\7"),
        DFA.unpack(u"\1\11"),
        DFA.unpack(u"\1\10\3\uffff\1\10\1\uffff\1\4\1\5\22\uffff\1\3\3"
        u"\uffff\1\7"),
        DFA.unpack(u"\1\10\3\uffff\1\10\31\uffff\1\7"),
        DFA.unpack(u"\1\10\3\uffff\1\10\1\uffff\1\11\1\12\26\uffff\1\7"),
        DFA.unpack(u""),
        DFA.unpack(u""),
        DFA.unpack(u"\1\10\3\uffff\1\10\1\uffff\1\11\1\12\26\uffff\1\7"),
        DFA.unpack(u"\1\10\3\uffff\1\10\31\uffff\1\7")
    ]

    # class definition for DFA #30

    DFA30 = DFA


    FOLLOW_cif_in_parse50 = frozenset([])
    FOLLOW_EOF_in_parse52 = frozenset([1])
    FOLLOW_COMMENTS_in_cif65 = frozenset([1, 5, 9])
    FOLLOW_WHITESPACE_in_cif70 = frozenset([1, 5, 9])
    FOLLOW_data_block_in_cif76 = frozenset([1, 5, 9])
    FOLLOW_WHITESPACE_in_cif80 = frozenset([5, 9])
    FOLLOW_data_block_in_cif83 = frozenset([1, 5, 9])
    FOLLOW_WHITESPACE_in_cif89 = frozenset([1])
    FOLLOW_value_in_loop_body112 = frozenset([1, 5])
    FOLLOW_WHITESPACE_in_loop_body125 = frozenset([5, 11, 13, 14, 31, 32, 33, 34])
    FOLLOW_value_in_loop_body139 = frozenset([1, 5])
    FOLLOW_SAVE_FRAME_HEADING_in_save_frame174 = frozenset([5])
    FOLLOW_WHITESPACE_in_save_frame178 = frozenset([5, 8, 10])
    FOLLOW_data_items_in_save_frame181 = frozenset([5])
    FOLLOW_WHITESPACE_in_save_frame186 = frozenset([5, 7])
    FOLLOW_SAVE_in_save_frame189 = frozenset([1])
    FOLLOW_TAG_in_data_items213 = frozenset([5])
    FOLLOW_WHITESPACE_in_data_items215 = frozenset([11, 13, 14, 31, 32, 33, 34])
    FOLLOW_value_in_data_items217 = frozenset([1])
    FOLLOW_loop_header_in_data_items236 = frozenset([5, 11, 13, 14, 31, 32, 33, 34])
    FOLLOW_WHITESPACE_in_data_items238 = frozenset([5, 11, 13, 14, 31, 32, 33, 34])
    FOLLOW_loop_body_in_data_items241 = frozenset([1])
    FOLLOW_DATA_BLOCK_HEADING_in_data_block258 = frozenset([1, 5])
    FOLLOW_WHITESPACE_in_data_block271 = frozenset([5, 6, 8, 10])
    FOLLOW_data_items_in_data_block276 = frozenset([1, 5])
    FOLLOW_save_frame_in_data_block280 = frozenset([1, 5])
    FOLLOW_LOOP__in_loop_header296 = frozenset([5])
    FOLLOW_WHITESPACE_in_loop_header304 = frozenset([5, 8])
    FOLLOW_TAG_in_loop_header307 = frozenset([5])
    FOLLOW_WHITESPACE_in_loop_header312 = frozenset([1])
    FOLLOW_31_in_inapplicable327 = frozenset([1])
    FOLLOW_32_in_unknown336 = frozenset([1])
    FOLLOW_inapplicable_in_value358 = frozenset([1])
    FOLLOW_unknown_in_value362 = frozenset([1])
    FOLLOW_33_in_value366 = frozenset([1])
    FOLLOW_char_string_in_value370 = frozenset([1])
    FOLLOW_numeric_in_value375 = frozenset([1])
    FOLLOW_text_field_in_value378 = frozenset([1])
    FOLLOW_DIGIT_in_unsigned_integer399 = frozenset([1, 11])
    FOLLOW_set_in_integer411 = frozenset([11, 33, 34])
    FOLLOW_unsigned_integer_in_integer422 = frozenset([1])
    FOLLOW_integer_in_float_432 = frozenset([12])
    FOLLOW_EXPONENT_in_float_434 = frozenset([1])
    FOLLOW_set_in_float_440 = frozenset([11, 31])
    FOLLOW_DIGIT_in_float_454 = frozenset([11, 31])
    FOLLOW_31_in_float_458 = frozenset([11, 33, 34])
    FOLLOW_unsigned_integer_in_float_460 = frozenset([1, 12])
    FOLLOW_DIGIT_in_float_466 = frozenset([11, 31])
    FOLLOW_31_in_float_470 = frozenset([1, 12])
    FOLLOW_EXPONENT_in_float_475 = frozenset([1])
    FOLLOW_integer_in_number487 = frozenset([1])
    FOLLOW_float__in_number491 = frozenset([1])
    FOLLOW_number_in_numeric500 = frozenset([1])
    FOLLOW_number_in_numeric506 = frozenset([35])
    FOLLOW_35_in_numeric508 = frozenset([11])
    FOLLOW_DIGIT_in_numeric511 = frozenset([11, 36])
    FOLLOW_36_in_numeric515 = frozenset([1])
    FOLLOW_CHAR_STRING_in_char_string527 = frozenset([1])
    FOLLOW_SEMI_COLON_TEXT_FIELD_in_text_field537 = frozenset([1])



def main(argv, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr):
    from antlr3.main import ParserMain
    main = ParserMain("cifProtoLexer", cifProtoParser)
    main.stdin = stdin
    main.stdout = stdout
    main.stderr = stderr
    main.execute(argv)


if __name__ == '__main__':
    main(sys.argv)
