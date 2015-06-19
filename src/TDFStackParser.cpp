#include <fstream>
#include <iostream>
#include <map>
#include <list>

#include "TDFStackParser.hpp"

// parse modes
const int M_NONE				= 0;
const int M_PARSE_SECTION_NAME	= 1, M_PARSE_SECTION		= 2;
const int M_PARSE_SECTION_KEY	= 3, M_PARSE_SECTION_VAL	= 4;
const int M_COMMENT_SL			= 5, M_COMMENT_ML			= 6;

const std::string parseModes[] = {
	"NONE",
	"PARSE_SECTION_NAME", "PARSE_SECTION",
	"PARSE_SECTION_KEY", "PARSE_SECTION_VAL"
	"COMMENT_SL", "COMMENT_ML",
};

const char sep = '\\';
const std::string tokens = "{}[];=/*\n";

static void ReplaceStringChars(std::string& s, char c1, char c2) {
	for (unsigned int i = 0; i < s.size(); i++) {
		if (s[i] == c1) {
			s[i] = c2;
		}
	}
}



CTDFStackParser::CTDFStackParser(const std::string& fname) {
	if (LoadFile(fname)) {
		Parse();
	}
}



void CTDFStackParser::PrintModeStack() const {
	for (std::list<int>::const_iterator it = modeStack.begin(); it != modeStack.end(); it++) {
		std::cout << parseModes[*it] << " ";
	}

	std::cout << std::endl;
}

bool CTDFStackParser::LoadFile(const std::string& fname) {
	std::ifstream stream(fname.c_str(), std::ios::in);
	std::string line;

	if (!stream.good()) {
		return false;
	}

	while (!stream.eof()) {
		std::getline(stream, line);
		fileString += (line + '\n');
	}

	stream.close();
	return true;
}

bool CTDFStackParser::Parse() {
	ReplaceStringChars(fileString, '\r', '\n');
	ReplaceStringChars(fileString, '\t', ' ');

	std::string sectionStackKey;

	std::string sec; // section name currently being read
	std::string key; // key currently being read
	std::string val; // value currently being read

	sectionStack.clear();
	modeStack.push_back(M_NONE);

	unsigned int k = fileString.size();

	for (unsigned int i = 0; i < k; i++) {
		unsigned int j = i + 1;
		unsigned char c = tolower(fileString[i]);
		unsigned int t = tokens.find(c);

		if (t != std::string::npos) {
			const int m = modeStack.back();

			if (c == '[') {
				// '[' can occur in values too
				if (m == M_PARSE_SECTION_VAL) {
					val += c;
				} else {
					if (m == M_NONE || m == M_PARSE_SECTION) {
						modeStack.push_back(M_PARSE_SECTION_NAME);
					} else {
						parseLog += "syntax error: stray '[' encountered\n";
					}
				}
			}
			if (c == ']') {
				// ']' can occur in values too
				if (m == M_PARSE_SECTION_VAL) {
					val += c;
				} else {
					if (m == M_PARSE_SECTION_NAME) {
						if (sectionStack.empty())
							sectionStackKey += sec;
						else
							sectionStackKey += (sep + sec);

						modeStack.pop_back();
						sectionStack.push_back(sec);
						sec.clear();
					} else {
						parseLog += "syntax error: stray ']' encountered\n";
					}
				}
			}

			if (c == '{') {
				// '{' can occur in values too
				if (m == M_PARSE_SECTION_VAL) {
					val += c;
				} else {
					if (m == M_NONE || m == M_PARSE_SECTION) {
						// we want to read a section name or key now
						modeStack.push_back(M_PARSE_SECTION);
					} else {
						parseLog += "syntax error: stray '{' encountered\n";
					}
				}
			}
			if (c == '}') {
				// '}' can occur in values too
				if (m == M_PARSE_SECTION_VAL) {
					val += c;
				} else {
					if (m == M_PARSE_SECTION) {
						sectionStack.pop_back();
						modeStack.pop_back();

						j = sectionStackKey.find_last_of(sep);
						sectionStackKey = sectionStackKey.substr(0, j);
					} else {
						parseLog += "syntax error: stray '}' encountered\n";
					}
				}
			}

			if (c == '=') {
				if (m == M_PARSE_SECTION_KEY) {
					// we want to read a section value now
					modeStack.pop_back();
					modeStack.push_back(M_PARSE_SECTION_VAL);
				} else {
					if (m != M_COMMENT_SL && m != M_COMMENT_ML) {
						parseLog += "syntax error: stray '=' encountered\n";
					}
				}
			}
			if (c == ';') {
				if (m == M_PARSE_SECTION_VAL) {
					modeStack.pop_back();

					keyVals[sectionStackKey][key] = val;
					// keyVals[sectionStack.back()][key] = val;
					key.clear();
					val.clear();
				} else {
					if (m != M_COMMENT_SL && m != M_COMMENT_ML) {
						parseLog += "syntax error: stray ';' encountered\n";
					}
				}
			}


			if (c == '/') {
				// '/' can occur in values too
				if (m == M_PARSE_SECTION_VAL) {
					val += c;
				} else {
					if (i < k - 1) {
						if (fileString[j] == '/') {
							if (m != M_COMMENT_SL && m != M_COMMENT_ML) {
								modeStack.push_back(M_COMMENT_SL);
							}
						}
						if (fileString[j] == '*') {
							if (m != M_COMMENT_SL && m != M_COMMENT_ML) {
								modeStack.push_back(M_COMMENT_ML);
							}
						}
					}
				}
			}

			if (c == '*') {
				// '*' can occur in values too
				if (m == M_PARSE_SECTION_VAL) {
					val += c;
				} else {
					if (i < k - 1 && fileString[j] == '/') {
						if (m == M_COMMENT_ML) {
							modeStack.pop_back();
						} else {
							parseLog += "syntax error: stray */ comment terminator encountered\n";
						}
					}
				}
			}

			if (c == '\n') {
				if (m == M_COMMENT_SL) {
					modeStack.pop_back();
				}
			}
		} else {
			// character is not a token
			if (modeStack.back() == M_PARSE_SECTION && c != ' ') {
				modeStack.push_back(M_PARSE_SECTION_KEY);
			}

			// FIXME: only section values may contain spaces
			if (modeStack.back() == M_PARSE_SECTION_NAME) { sec += c; }
			if (modeStack.back() == M_PARSE_SECTION_KEY) { key += c; }
			if (modeStack.back() == M_PARSE_SECTION_VAL) { val += c; }
		}
	}

	return (modeStack.size() == 1 && modeStack.back() == M_NONE && parseLog.size() == 0);
}

std::string CTDFStackParser::GetParseResult() {
	std::string s;

	for (SectionKeyValMap::iterator it = keyVals.begin(); it != keyVals.end(); it++) {
		const std::string& section = it->first;
		const std::map<std::string, std::string>& sectionKVPairs = keyVals[section];

		s += (section + '\n');

		for (std::map<std::string, std::string>::const_iterator p = sectionKVPairs.begin(); p != sectionKVPairs.end(); p++) {
			s += ('\t' + p->first + '=' + p->second + '\n');
		}
	}

	return s;
}
