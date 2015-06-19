#ifndef TDFSTACKPARSER_HDR
#define TDFSTACKPARSER_HDR

#include <iostream>
#include <sstream>
#include <map>

#include "vec3.hpp"

// simple stack-based TDF parser
class CTDFStackParser {
	public:
		CTDFStackParser(const std::string& = "");
		~CTDFStackParser() {}

		bool LoadFile(const std::string&);
		bool Parse();

		std::string GetValueForKey(const std::string& section, const std::string& key, const std::string& def) {
			if (keyVals.find(section) != keyVals.end()) {
				if (keyVals[section].find(key) != keyVals[section].end()) {
					return keyVals[section][key];
				}
			}

			return def;
		}

		// NOTE: not templated (would need to overload operator >>)
		vec3 GetVec(const std::string& section, const std::string& key) {
			vec3 v;

			std::stringstream ss;
			ss << GetValueForKey(section, key, "");
			ss >> v.x;
			ss >> v.y;
			ss >> v.z;

			return v;
		}
		template<typename T> T GetValue(const std::string& section, const std::string& key) {
			T t;

			std::stringstream ss;
			ss << GetValueForKey(section, key, "0");
			ss >> t;

			return t;
		}

		const std::string& GetParseLog() const { return parseLog; }
		std::string GetParseResult();

	private:
		void PrintModeStack() const;

		std::string fileString;
		std::string parseLog;

		// {"section name 1": {"key 1": "value 1", "key 2": "value 2" ...}}
		typedef std::map<std::string, std::map<std::string, std::string> > SectionKeyValMap;
		SectionKeyValMap keyVals;

		std::list<std::string> sectionStack;
		std::list<int> modeStack;
};

#endif
