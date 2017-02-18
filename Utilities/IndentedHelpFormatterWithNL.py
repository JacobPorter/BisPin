from optparse import IndentedHelpFormatter
import textwrap

""""
Code taken from https://groups.google.com/forum/#!topic/comp.lang.python/bfbmtUGhW8I
Also see https://groups.google.com/forum/#!topic/comp.lang.python/5y3u53nZmJs
Some discussion at https://bytes.com/topic/python/answers/734066-how-output-newline-carriage-return-optparse
@author: Dan Williams
@author: Tim Chase
"""

class IndentedHelpFormatterWithNL(IndentedHelpFormatter):
    """
    This formatter for optparse is able to add new lines into help text for an option.
    Example: p = optparse.OptionParser(formatter=IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    """
    def format_description(self, description):
        if not description:
            return ""
        desc_width = self.width - self.current_indent
        indent = " " * self.current_indent
        bits = description.split('\n')
        formatted_bits = [textwrap.fill(bit, desc_width, initial_indent=indent, subsequent_indent=indent) for bit in bits]
        result = "\n".join(formatted_bits) + "\n"
        return result

    def format_option(self, option):
        """
        The help for each option consists of two parts:
          * the opt strings and metavars
          eg. ("-x", or "-fFILENAME, --file=FILENAME")
          * the user-supplied help string
          eg. ("turn on expert mode", "read data from FILENAME")
        
        If possible, we write both of these on the same line:
          -x    turn on expert mode
        
        But if the opt string list is too long, we put the help
        string on a second line, indented to the same column it would
        start in if it fit on the first line.
          -fFILENAME, --file=FILENAME
              read data from FILENAME
          """
        result = []
        opts = self.option_strings[option]
        opt_width = self.help_position - self.current_indent - 2
        if len(opts) > opt_width:
            opts = "%*s%s\n" % (self.current_indent, "", opts)
            indent_first = self.help_position
        else: # start help on same line as opts
            opts = "%*s%-*s  " % (self.current_indent, "", opt_width, opts)
            indent_first = 0
        result.append(opts)
        if option.help:
            #help_text = option.help
            help_text = self.expand_default(option)
        # Everything is the same up through here
            help_lines = []
            for para in help_text.split("\n"):
                help_lines.extend(textwrap.wrap(para, self.help_width))
#             help_text = "\n".join([x.strip() for x in help_text.split("\n")])
#             for para in help_text.split("\n\n"):
#                 help_lines.extend(textwrap.wrap(para, self.help_width))
#                 if len(help_lines):
#                     # for each paragraph, keep the double newlines..
#                     help_lines[-1] += "\n"
            # Everything is the same after here
            result.append("%*s%s\n" % (indent_first, "", help_lines[0]))
            result.extend(["%*s%s\n" % (self.help_position, "", line) for line in help_lines[1:]])
        elif opts[-1] != "\n":
            result.append("\n")
        return "".join(result)