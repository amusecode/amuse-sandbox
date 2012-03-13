%% 
%% This is file `iopbk2e.cls'
%% 
%% Class file for Books to be published by Institute of Physics Publishing
%%
%% This file is distributed in the hope that it will be useful, 
%% but WITHOUT ANY WARRANTY; without even the implied warranty of 
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%% 
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{iopbk2e}[1997/2/13  IOP book style class]
\newcommand\@ptsize{0} 
\newif\if@restonecol
\newif\if@titlepage
\newif\ifiopams
\@titlepagefalse
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\DeclareOption{Crown}
   {\setlength\textwidth{35pc}
    \setlength\textheight{48pc}}
\DeclareOption{Demi}
   {\setlength\textwidth{24pc}
    \setlength\textheight{40pc}}
\DeclareOption{Royal}
   {\setlength\textwidth{27pc}
    \setlength\textheight{44pc}}
\DeclareOption{10pt}{\renewcommand\@ptsize{0}}
\DeclareOption{11pt}{\renewcommand\@ptsize{2}} % No 11pt version
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\DeclareOption{titlepage}{\@titlepagetrue}
\DeclareOption{notitlepage}{\@titlepagefalse}
\ExecuteOptions{a4paper,final,Royal}
\ProcessOptions
\DeclareMathAlphabet{\bi}{OML}{cmm}{b}{it}
\input{lib/iopbk1\@ptsize.clo}
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0\p@ \@plus \p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setlength\parindent{2em}
\setcounter{topnumber}{8}
\renewcommand\topfraction{1}
\setcounter{bottomnumber}{3}
\renewcommand\bottomfraction{.8}
\setcounter{totalnumber}{8}
\renewcommand\textfraction{0.1}
\renewcommand\floatpagefraction{.8}
\setcounter{dbltopnumber}{6}
\renewcommand\dbltopfraction{1}
\renewcommand\dblfloatpagefraction{.8}
%
%
% Headings for all pages apart from first
%
\def\ps@headings{\let\@oddfoot\@empty
      \let\@evenfoot\@empty
      \def\@evenhead{\thepage\hfil\itshape\rightmark}%
      \def\@oddhead{{\itshape\leftmark}\hfil\thepage}%
      \let\@mkboth\markboth
      \let\sectionmark\@gobble
      \let\subsectionmark\@gobble}
%
% Headings for first page
%
\def\ps@myheadings{\let\@oddfoot\@empty\let\@evenfoot\@empty
    \let\@oddhead\@empty\let\@evenhead\@empty
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble}
%
\def\maketitle{\begin{titlepage}
 \let\footnotesize\small \let\footnoterule\relax \setcounter{page}{0}
 \null
 \vskip 30pt 
 \begin{center}
 {\LARGE\bf \@title \@@par} \vskip 3em {\lineskip .75em
 \def\tabular{\Large\def\@halignto{}\@tabular}
 \begin{tabular}[t]{c}\@author 
 \end{tabular}\@@par} 
 \vskip 1.5em {\large \@date \@@par} \end{center} \@@par
\@thanks
\vfil                
\null
\end{titlepage}
\setcounter{footnote}{0} \let\thanks\relax
\gdef\@thanks{}\gdef\@author{}\gdef\@title{}\let\maketitle\relax} 

\def\labelenumi{(\roman{enumi})} 
\def\theenumi{\roman{enumi}} 
\def\labelenumii{(\alph{enumii})}
\def\theenumii{\alph{enumii}}
\def\p@enumii{\theenumi.}
\def\labelenumiii{\arabic{enumiii}.}
\def\theenumiii{\arabic{enumiii}}
\def\p@enumiii{\p@enumii.\theenumii}
\def\labelenumiv{\arabic{enumiv}.}
\def\theenumiv{\arabic{enumiv}} 
\def\p@enumiv{\p@enumiii.\theenumiii}

\def\labelitemi{$\bullet$}
\def\labelitemii{\bf --}
\def\labelitemiii{$\ast$}
\def\labelitemiv{$\cdot$}

\def\verse{\let\\=\@centercr 
 \list{}{\itemsep\z@ \itemindent -1.5em\listparindent \itemindent 
 \rightmargin\leftmargin\advance\leftmargin 1.5em}\item[]}
\let\endverse\endlist
\def\quotation{\list{}{\listparindent 1.5em
 \itemindent\listparindent
 \rightmargin\leftmargin\parsep 0pt plus 1pt}\item[]}
\let\endquotation=\endlist
\def\quote{\list{}{\rightmargin\leftmargin}\item[]}
\let\endquote=\endlist

\def\descriptionlabel#1{\hspace\labelsep \bf #1}
\def\description{\list{}{\labelwidth\z@ \itemindent-\leftmargin
 \let\makelabel\descriptionlabel}}
\let\enddescription\endlist
\newdimen\descriptionmargin
\descriptionmargin=3em

\def\enumerate{\ifnum \@enumdepth >3 \@toodeep\else
      \advance\@enumdepth \@ne 
      \edef\@enumctr{enum\romannumeral\the\@enumdepth}\list
      {\csname label\@enumctr\endcsname}{\usecounter
        {\@enumctr}\def\makelabel##1{##1\hss}%
    \labelsep=3\p@\parsep=0pt\listparindent=0pt\itemsep=0pt\partopsep=0pt
    \rightmargin=0pt}\fi}

\def\itemize{\ifnum \@itemdepth >3 \@toodeep\else \advance\@itemdepth \@ne
    \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
    \list{\csname\@itemitem\endcsname}{\def\makelabel##1{##1\hss}%
    \labelsep=3\p@\parsep=0pt\listparindent=0pt\itemsep=0pt\partopsep=0pt
    \rightmargin=0pt}\fi}

\def\titlepage{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
 \else \newpage \fi \thispagestyle{empty}\c@page\z@}

\def\endtitlepage{\if@restonecol\twocolumn \else \newpage \fi}

\arraycolsep 5pt 
\tabcolsep 6pt 
\extracolsep{0pt plus4pt minus 2pt}
\arrayrulewidth .4pt 
\doublerulesep 2pt 
\tabbingsep \labelsep 
\skip\@mpfootins = \skip\footins
\fboxsep = 3pt 
\fboxrule = .4pt 

\newcounter{part}
\newcounter {chapter}
\newcounter {section}[chapter]
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]

\def\thepart {\arabic{part}}
\def\thechapter {\arabic{chapter}}
\def\thesection {\thechapter.\arabic{section}}
\def\thesubsection {\thesection.\arabic{subsection}}
\def\thesubsubsection {\thesubsection .\arabic{subsubsection}}
\def\theparagraph {\thesubsubsection.\arabic{paragraph}}
\def\thesubparagraph {\theparagraph.\arabic{subparagraph}}
\def\@chapapp{Chapter}


\def\@pnumwidth{1.55em}
\def\@tocrmarg {2.55em}
\def\@dotsep{4.5}
\setcounter{tocdepth}{2}


\def\tableofcontents{\@restonecolfalse\if@twocolumn\@restonecoltrue
 \onecolumn\fi\chapter*{Contents}{}\thispagestyle{empty}
 \@starttoc{toc}\if@restonecol\twocolumn\fi}

\def\l@part#1#2{\addpenalty{-\@highpenalty}
 \addvspace{2.25em plus 1pt} \begingroup
 \@tempdima 1.5em \parindent \z@ \rightskip \@pnumwidth \parfillskip
-\@pnumwidth 
 {\bf \leavevmode\leftskip=\@tempdima #1\hfil 
  \hbox to\@pnumwidth{\hss #2}\par}
 \nobreak \endgroup}

\def\l@chapter#1#2{\pagebreak[3] 
 \vskip 0.5em plus 1pt \@tempdima 1.5em \begingroup
 \parindent \z@ \rightskip \@pnumwidth 
 \parfillskip -\@pnumwidth 
 \bf \leavevmode #1\hfil \hbox to\@pnumwidth{\hss #2}\par
 \endgroup}
%
\def\l@unnumbered#1#2{\pagebreak[3] 
 \vskip 0.5em plus 1pt \@tempdima 1.5em \begingroup
 \parindent \z@ \rightskip \@pnumwidth 
 \parfillskip -\@pnumwidth 
 \bf \leavevmode #1\hfil \hbox to\@pnumwidth{\hss #2}\par
 \endgroup}
%
\def\l@section{\@dottedtocline{1}{1.5em}{2.3em}}
\def\l@subsection{\@dottedtocline{2}{3.8em}{3.2em}}
\def\l@subsubsection{\@dottedtocline{3}{7.0em}{4.1em}}
\def\l@paragraph{\@dottedtocline{4}{10em}{5em}}
\def\l@subparagraph{\@dottedtocline{5}{12em}{6em}}
\def\listoffigures{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
 \fi\chapter*{List of Figures\@mkboth
 {LIST OF FIGURES}{LIST OF FIGURES}}\@starttoc{lof}\if@restonecol\twocolumn
 \fi}
\def\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\def\listoftables{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
 \fi\chapter*{List of Tables\@mkboth
 {LIST OF TABLES}{LIST OF TABLES}}\@starttoc{lot}\if@restonecol\twocolumn
 \fi}
\let\l@table\l@figure
%
% Redefinition to remove dotted lines from \@dottedtocline
%
\def\@dottedtocline#1#2#3#4#5{\ifnum #1>\c@tocdepth \else
  \vskip \z@ plus .2pt
  {\leftskip #2\relax \rightskip \@tocrmarg \parfillskip -\rightskip
    \parindent #2\relax\@afterindenttrue
   \interlinepenalty\@M
   \leavevmode
   \@tempdima #3\relax \advance\leftskip \@tempdima \hbox{}\hskip -\leftskip
    #4\nobreak\hfill \nobreak \hbox to\@pnumwidth{\hfil\rm #5}\par}\fi}
 

\newcommand{\Bibliography}{\chapter*{Bibliography}}
\newcommand{\References}{\chapter*{References}}

\def\thebibliography#1{\list
 {\hfil[\arabic{enumi}]}{\labelsep=5pt\itemindent=-10pt
 \parsep=0pt\itemsep=0pt\settowidth\labelwidth{\small[#1]}%
 \leftmargin\labelwidth
 \advance\leftmargin\labelsep
 \advance\leftmargin -\itemindent
 \usecounter{enumi}}\small
 \def\newblock{\hskip .11em plus .33em minus .07em}
 \sloppy\clubpenalty4000\widowpenalty4000
 \sfcode`\.=1000\relax}
\let\endthebibliography=\endlist
%
% Macro to used for references in the Harvard system
%
\newenvironment{thereferences}{\begin{thereferences}}{\end{thereferences}}
%
\def\thereferences{\list{}{\itemsep=0pt\labelsep=0pt\itemindent=-20pt
\parsep=0pt\labelwidth=0pt\leftmargin=20pt}\small
\def\newblock{\hskip .11em plus .33em minus .07em}
\sloppy\clubpenalty4000\widowpenalty4000
\sfcode`\.=1000\relax}
\let\endthereferences=\endlist
%
 

\newif\if@restonecol
% sll
%\def\theindex{\@restonecoltrue\if@twocolumn\@restonecolfalse\fi
%\columnseprule \z@
%\columnsep 35pt
%\parindent\z@
%\parskip\z@ plus .3pt\relax
%\let\item\@idxitem
%\twocolumn[\chapter*{Index}]}
%\def\@idxitem{\par\hangindent 30pt}
%\def\subitem{\par\hangindent 30pt \hspace*{10pt}}
%\def\subsubitem{\par\hangindent 30pt \hspace*{20pt}}
%\def\endtheindex{\if@restonecol\onecolumn\else\clearpage\fi}
%\def\indexspace{\par \vskip 10pt plus 5pt minus 3pt\relax}
\newcommand{\indexname}{Index}
\newenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \columnseprule \z@
                \columnsep 35\p@
                \twocolumn[\@makeschapterhead{\indexname}]%
%                \twocolumn[\@makeschapterhead{Index}]%
%                \@mkboth{\uppercase{\indexname}}%
%                        {\uppercase{\indexname}}%
                \@mkboth{\indexname}%
                        {\indexname}%
                \parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand{\@idxitem}{\par\hangindent 30pt}
\newcommand{\subitem}{\par\hangindent 30pt \hspace*{10pt}}
\newcommand{\subsubitem}{\par\hangindent 30pt \hspace*{20pt}}
\newcommand{\indexspace}{\par \vskip 10pt plus 5pt minus 3pt\relax}
% sll

\def\footnoterule{}%
\@addtoreset{footnote}{chapter} 
\long\def\@makefntext#1{\parindent 1em\noindent 
 \makebox[1em][l]{$^{\@thefnmark}$}#1}


\setcounter{topnumber}{3}
\def\topfraction{1}
\setcounter{bottomnumber}{2}
\def\bottomfraction{.8}
\setcounter{totalnumber}{3}
\def\textfraction{0}
\def\floatpagefraction{.8}
\setcounter{dbltopnumber}{3}
\def\dbltopfraction{1}
\def\dblfloatpagefraction{.8}

\long\def\@makecaption#1#2{\vskip 10pt 
 \setbox\@tempboxa\hbox{\small{\bf #1.} #2}%
 \ifdim \wd\@tempboxa >\hsize \unhbox\@tempboxa\par \else \hbox
to\hsize{\hfil\box\@tempboxa\hfil} 
 \fi}

\def\tabular{\small\def\@halignto{}\@tabular}

\newcounter{figure}[chapter]
\def\thefigure{\thechapter.\@arabic\c@figure}
\def\fps@figure{htbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{Figure \thefigure}
\def\figure{\@float{figure}}
\let\endfigure\end@float
\@namedef{figure*}{\@dblfloat{figure}}
\@namedef{endf*}{\end@dblfloat}
\newcounter{table}[chapter]
\def\thetable{\thechapter.\@arabic\c@table}
\def\fps@table{htbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{Table \thetable}
\def\table{\@float{table}}
\let\endtable\end@float
\@namedef{table*}{\@dblfloat{table}}
\@namedef{endtable*}{\end@dblfloat}

\def\maketitle{\begin{titlepage}
 \let\footnotesize\small \let\footnoterule\relax \setcounter{page}{0}
 \null
 \vskip 30pt 
 \begin{center}
 \Large
 {\LARGE\bf \@title \par} \vskip 3em {\lineskip .75em
 \def\tabular{\Large\def\@halignto{}\@tabular}
 \begin{tabular}[t]{c}\@author 
 \end{tabular}\par} 
 \vskip 1.5em {\large \@date \par} \end{center} \par
\@thanks
\vfil                
\null
\end{titlepage}
\setcounter{footnote}{0} \let\thanks\relax
\gdef\@thanks{}\gdef\@author{}\gdef\@title{}\let\maketitle\relax} 


\mark{{}{}} 

\def\ps@headings{\let\@mkboth\markboth
 \def\@oddfoot{}%
 \def\@evenfoot{}%
 \def\@evenhead{\makebox[4pc][l]{\normalsize\rm \thepage}%
  \normalsize\sl\leftmark\hfill}%
 \def\@oddhead{\hfill{\normalsize\sl\rightmark}%
  \makebox[4pc][r]{\normalsize\rm\thepage}}%
 \def\chaptermark##1{\markboth{##1}{##1}}%
 \def\sectionmark##1{\markright{##1}}}%

\def\ps@myheadings{\let\@mkboth\markboth
 \def\@oddhead{}%
 \def\@oddfoot{\hfill\normalsize\rm\thepage}%
 \def\@evenhead{}%
 \def\@evenfoot{\normalsize\rm\thepage\hfill}%
 \def\chaptermark##1{\markboth{##1}{##1}}%
 \def\sectionmark##1{\markright{##1}}}


\def\today{\ifcase\month\or
 January\or February\or March\or April\or May\or June\or
 July\or August\or September\or October\or November\or December\fi
 \space\number\day, \number\year}
\@addtoreset{equation}{chapter} 
\def\theequation{\thechapter.\arabic{equation}}

\def\@begintheorem#1#2{\rm \trivlist \item[\hskip \labelsep{\it #1\ #2.}]}
\def\@opargbegintheorem#1#2#3{\rm \trivlist
      \item[\hskip \labelsep{\it #1\ #2\ (#3).}]}
%
\def\eqnarray{\stepcounter{equation}\let\@currentlabel=\theequation
\global\@eqnswtrue
\global\@eqcnt\z@\tabskip\@centering\let\\=\@eqncr
$$\halign to \displaywidth\bgroup\@eqnsel\hskip\@centering
  $\displaystyle\tabskip\z@{##{}}$&\global\@eqcnt\@ne 
  \hfil$\displaystyle{{}##{}}$\hfil
  &\global\@eqcnt\tw@ $\displaystyle\tabskip\z@{{}##}$\hfil 
   \tabskip\@centering&\llap{##}\tabskip\z@\cr}
%
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\ifiopams
\renewcommand{\cal}{\protect\pcal}
\else
\newcommand{\cal}{\protect\pcal}
\fi
\newcommand{\pcal}{\@fontswitch{\relax}{\mathcal}}
\ifiopams
\renewcommand{\mit}{\protect\pmit}
\else
\newcommand{\mit}{\protect\pmit}
\fi
\newcommand{\pmit}{\@fontswitch{\relax}{\mathnormal}}
% sll 
\newcommand{\text}[1]{\mbox{#1}}
\newcommand{\bold}[1]{{\bf #1}}
%
\newcommand{\e}{\mathrm{e}}
\newcommand{\rme}{\mathrm{e}}
\newcommand{\rmi}{\mathrm{i}}
\newcommand{\rmd}{\mathrm{d}}
\newcommand{\Tr}{\mathop{\mathrm{Tr}}\nolimits}
\newcommand{\tr}{\mathop{\mathrm{tr}}\nolimits}
\newcommand{\Or}{\mathop{\mathrm{Or}}\nolimits}
\newcommand{\tdot}[1]{\stackrel{\dots}{#1}}   % Added 1/9/94
\newcommand{\lshad}{[\![}
\newcommand{\rshad}{]\!]}
\newcommand{\case}[2]{{\textstyle\frac{#1}{#2}}}
\def\pt(#1){({\it #1\/})}
\newcommand{\dsty}{\displaystyle}
\newcommand{\tsty}{\textstyle}
\newcommand{\ssty}{\scriptstyle}
\newcommand{\sssty}{\scriptscriptstyle}
%
\def\;{\protect\psemicolon}
\def\psemicolon{\relax\ifmmode\mskip\thickmuskip\else\kern .3333em\fi}
%
\newcommand{\eref}[1]{(\ref{#1})}
\newcommand{\sref}[1]{section~\ref{#1}}
\newcommand{\cref}[1]{chapter~\ref{#1}}
\newcommand{\fref}[1]{figure~\ref{#1}}
\newcommand{\tref}[1]{table~\ref{#1}}
\newcommand{\Eref}[1]{Equation (\ref{#1})}
\newcommand{\Sref}[1]{Section~\ref{#1}}
\newcommand{\Cref}[1]{Chapter~\ref{#1}}
\newcommand{\Fref}[1]{Figure~\ref{#1}}
\newcommand{\Tref}[1]{Table~\ref{#1}}

\newcommand{\opencircle}{\mbox{\Large$\circ\,$}}  % moved Large outside maths
\newcommand{\opensquare}{\mbox{$\rlap{$\sqcap$}\sqcup$}}
\newcommand{\opentriangle}{\mbox{$\triangle$}}
\newcommand{\opentriangledown}{\mbox{$\bigtriangledown$}}
\newcommand{\opendiamond}{\mbox{$\diamondsuit$}}
\newcommand{\fullcircle}{\mbox{{\Large$\bullet\,$}}} % moved Large outside maths
\newcommand{\fullsquare}{\,\vrule height5pt depth0pt width5pt}
\newcommand{\dotted}{\mbox{${\mathinner{\cdotp\cdotp\cdotp\cdotp\cdotp\cdotp}}$}}
\newcommand{\dashed}{\mbox{-\; -\; -\; -}}
\newcommand{\broken}{\mbox{-- -- --}}
\newcommand{\longbroken}{\mbox{--- --- ---}}
\newcommand{\chain}{\mbox{--- $\cdot$ ---}}
\newcommand{\dashddot}{\mbox{--- $\cdot$ $\cdot$ ---}}
\newcommand{\full}{\mbox{------}}
%
\newcommand{\figblank}{\begin{figure}[p]
 \vspace*{44pc}
 \addtocounter{figure}{1}
 \end{figure}}
%
\newcommand{\tabblank}{\begin{table}[p]
 \vspace*{44pc}
 \addtocounter{table}{1}
 \end{table}}

%%%%%%%%%%%%%%%%%%%%%
% Tables rules      %
%%%%%%%%%%%%%%%%%%%%%

\newcommand{\boldarrayrulewidth}{1pt} 
% Width of bold rule in tabular environment.

\def\bhline{\noalign{\ifnum0=`}\fi\hrule \@height  
\boldarrayrulewidth \futurelet \@tempa\@xhline}

\def\@xhline{\ifx\@tempa\hline\vskip \doublerulesep\fi
      \ifnum0=`{\fi}}

%
% Rules for tables with extra space around
%
\newcommand{\br}{\ms\bhline\ms}
\newcommand{\mr}{\ms\hline\ms}
\newcommand{\ms}{\noalign{\vspace{3pt plus2pt minus1pt}}}
\newcommand{\bs}{\noalign{\vspace{6pt plus2pt minus2pt}}}
\newcommand{\ns}{\noalign{\vspace{-3pt plus-1pt minus-1pt}}}
\newcommand{\centre}[2]{\multicolumn{#1}{c}{#2}}
\newcommand{\crule}[1]{\multispan{#1}{\hrulefill}}
%
% Macros for aligning columns
%
\newcommand{\0}{\phantom{0}}
\newcommand{\m}{\llap{$-$}}
\newcommand{\phm}{\phantom{$-$}}
%
\newsavebox{\tablebox}
\newcommand{\notedtable}[1]{\savebox{\tablebox}%
{#1}\usebox{\tablebox}}
\newcommand{\tabnote}[1]{\parbox{\wd\tablebox}{\strut\small #1}\par}
%
\def\verbatim{\@verbatim \small\frenchspacing\@vobeyspaces \@xverbatim}
\@namedef{verbatim*}{\@verbatim\small\@sxverbatim}
\def\endverbatim{\endtrivlist\normalsize\rm\relax}
\expandafter\let\csname endverbatim*\endcsname =\endverbatim
%
\newcommand{\etal}{{\em et al\/}\ }
\newcommand{\dash}{------}
\newcommand{\nonum}{\item[]}
%
% For index
%
\newcommand{\1}[1]{\item\uppercase{#1}}
\newcommand{\2}[1]{\subitem\uppercase{#1}}
%
\def\LaTeX{L\kern-.26em \raise.6ex\hbox{\tiny A}%
   \kern-.1em\TeX}%
\def\bigLaTeX{L\kern-.26em \raise.6ex\hbox{\small\bf A}%
   \kern-.1em\TeX}%
%
% abbreviations for IOPP journals
%
\newcommand{\CQG}{{\em Class. Quantum Grav.} }
\newcommand{\IP}{{\em Inverse Problems\/} }
\newcommand{\JPA}{{\em J. Phys. A: Math. Gen.} }
\newcommand{\JPB}{{\em J. Phys. B: At. Mol. Phys.} }      %1968-87
\newcommand{\jpb}{{\em J. Phys. B: At. Mol. Opt. Phys.} } %1988 and onwards
\newcommand{\JPC}{{\em J. Phys. C: Solid State Phys.} }   %1968--1988
\newcommand{\JPCM}{{\em J. Phys.: Condens. Matter\/} }    %1989 and onwards
\newcommand{\JPD}{{\em J. Phys. D: Appl. Phys.} }
\newcommand{\JPE}{{\em J. Phys. E: Sci. Instrum.} }
\newcommand{\JPF}{{\em J. Phys. F: Met. Phys.} }
\newcommand{\JPG}{{\em J. Phys. G: Nucl. Phys.} }         %1975--1988
\newcommand{\jpg}{{\em J. Phys. G: Nucl. Part. Phys.} }   %1989 and onwards
\newcommand{\NL}{{\em Nonlinearity\/} }
\newcommand{\PMB}{{\em Phys. Med. Biol.} }
\newcommand{\RPP}{{\em Rep. Prog. Phys.} }
\newcommand{\SST}{{\em Semicond. Sci. Technol.} }
\newcommand{\SUST}{{\em Supercond. Sci. Technol.} }
\newcommand{\MSMSE}{{\em Modelling Simulation Mater. Sci. Eng.} }
\newcommand{\MST}{{\em Meas. Sci. Technol.} }             %1990 and onwards
\newcommand{\NET}{{\em Network\/} }
\newcommand{\NT}{{\em Nanotechnology} }
\newcommand{\PAO}{{\em Pure Appl. Optics\/} }
\newcommand{\PSST}{{\em Plasma Sources Sci. Technol.} }
\newcommand{\QO}{{\em Quantum Opt.} }
\newcommand{\WRM}{{\em Waves Random Media\/} }
%
% Other commonly quoted journals
%
\newcommand{\AC}{{\em Acta Crystallogr.} }
\newcommand{\AM}{{\em Acta Metall.} }
\newcommand{\AP}{{\em Ann. Phys., Lpz.} }
\newcommand{\APNY}{{\em Ann. Phys., NY\/} }
\newcommand{\APP}{{\em Ann. Phys., Paris\/} }
\newcommand{\CJP}{{\em Can. J. Phys.} }
\newcommand{\JAP}{{\em J. Appl. Phys.} }
\newcommand{\JCP}{{\em J. Chem. Phys.} }
\newcommand{\JJAP}{{\em Japan. J. Appl. Phys.} }
\newcommand{\JP}{{\em J. Physique\/} }
\newcommand{\JPhCh}{{\em J. Phys. Chem.} }
\newcommand{\JMMM}{{\em J. Magn. Magn. Mater.} }
\newcommand{\JMP}{{\em J. Math. Phys.} }
\newcommand{\JOSA}{{\em J. Opt. Soc. Am.} }
\newcommand{\JPSJ}{{\em J. Phys. Soc. Japan\/} }
\newcommand{\JQSRT}{{\em J. Quant. Spectrosc. Radiat. Transfer\/} }
\newcommand{\NC}{{\em Nuovo Cimento\/} }
\newcommand{\NIM}{{\em Nucl. Instrum. Methods\/} }
\newcommand{\NP}{{\em Nucl. Phys.} }
\newcommand{\PL}{{\em Phys. Lett.} }
\newcommand{\PR}{{\em Phys. Rev.} }
\newcommand{\PRL}{{\em Phys. Rev. Lett.} }
\newcommand{\PRS}{{\em Proc. R. Soc.} }
\newcommand{\PS}{{\em Phys. Scr.} }
\newcommand{\PSS}{{\em Phys. Status Solidi\/} }
\newcommand{\PTRS}{{\em Phil. Trans. R. Soc.} }
\newcommand{\RMP}{{\em Rev. Mod. Phys.} }
\newcommand{\RSI}{{\em Rev. Sci. Instrum.} }
\newcommand{\SSC}{{\em Solid State Commun.} }
\newcommand{\ZP}{{\em Z. Phys.} }
%
\newif\ifnbk \nbkfalse
%
\ps@headings \pagenumbering{arabic} \onecolumn 
%
% Modified by Steve Lloyd 03/02/97 to fix index problems.
%
% Modified by Steve Lloyd 29/10/96 to change to a class file iopbk.cls.
%
%   - Add \NeedsTeXFormat \ProvidesClass
%   - Remove \typeout
%   - Add \newfam\bffam
%   - Change \@ptsize.sty to \@ptsize.clo
%   - Remove \newcommand{\mathrm}
%   - Add \DeclareOldFontCommand{\rm} etc
%
% sll 


