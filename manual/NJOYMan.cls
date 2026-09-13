% NJOY class file for typesetting the NJOY manual
% By Jeremy Lloyd Conlin jlconlin@lanl.gov
% October 26, 2011
%
% This class file is an attempt to upgrade the deprecated and extremely outdated
% file that is currently used.

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{NJOYMan}[2011/10/26 NJOY Manual class]

% Make this class a variant of the article class
\LoadClassWithOptions{article}

\RequirePackage{fancyvrb}
\RequirePackage{graphicx}
\RequirePackage[usenames, dvipsnames]{color}
\RequirePackage{hyperref}

\hypersetup{
    colorlinks=true,
    linkcolor=RoyalBlue,            % Make internal links blue
    citecolor=black,
    bookmarksopen=true,         % Open bookmarks by default
    bookmarksopenlevel=1        % but only to the first level
    bookmarksnumbered=true   % Include section numbers in bookmarks
    pdftitle={\@title},
    pdfauthor={\@author},
}

%-----------------------------------------------------------
% Define some commands to approximate the original NJOY
%\DefineShortVerb{\|}       % Verb environments only need to be delimited by '|' 
\newcommand{\cword}[1]{\Verb+#1+}

% Custom verbatim environment called ccode
\DefineVerbatimEnvironment{ccode}{Verbatim}{frame=single}

%-----------------------------------------------------------
% What follows extracted from Bob MacFarlane's la.sty file:
% (doesn't seem to hurt anything, but nothing seems different
%  either compared to not having this).
%    ****************************************
%    *              THE INDEX               *
%    ****************************************
%
% THE THEINDEX ENVIRONMENT
% Produces double column format, with each paragraph a separate entry.
% The user commands \item, \subitem and \subsubitem are used to
% produce the entries, and \indexspace adds an extra vertical space
% that's the right size to put above the first entry with a new letter
% of the alphabet.

\newif\if@restonecol

\def\theindex{\@restonecoltrue\if@twocolumn\@restonecolfalse\fi
\columnseprule \z@
\columnsep 35pt\twocolumn[\section*{Index}]
    \@mkboth{INDEX}{INDEX}\parindent\z@
    \parskip\z@ plus .3pt\relax\let\item\@idxitem
    \singlespacing}

\def\@idxitem{\par\hangindent 40pt}
\def\subitem{\par\hangindent 40pt \hspace*{20pt}}
\def\subsubitem{\par\hangindent 40pt \hspace*{30pt}}

\def\endtheindex{\if@restonecol\onecolumn\else\clearpage\fi}

\def\indexspace{\par \vskip 10pt plus 5pt minus 3pt\relax}
%-----------------------------------------------------------

