/*!
 * \file EventCut.h
 * \brief Base class for event cuts
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 1.0.1 $
 * \date $Date: 09/10/2024 $
 */

 #ifndef _EVENTCUT_H_
 #define _EVENTCUT_H_
 
 #include <iostream>
 #include <string>
 #include <vector>
 #include <utility>
 #include <cmath>  
 
 class PHCompositeNode;
 
 class EventCut
 {
   public:
 
     virtual ~EventCut() {}
     virtual void identify(std::ostream &os = std::cout) const { os << "EventCut::" << Name() << std::endl; }
     void clear(){ Passed(false); }
   
     // setters
     void Passed(bool result) { m_result = result; }
     void Name(const std::string &name) { m_name = name; }
 
     // getters
     bool Passed() const { return m_result; }
     const std::string& Name() const { return m_name; }
 
     int Verbosity() const { return m_verbosity; }
     void Verbosity(int verbosity) { m_verbosity = verbosity; }
 
     virtual bool operator()(PHCompositeNode* topNode) = 0;
 
     void SetNodeNames(const std::vector<std::string> & names) { m_node_names = names; }
     const std::vector<std::string>& GetNodeNames() const { return m_node_names; }
     void SetNodeName(const std::string &name) { SetNodeNames({ name }); }
     const std::string& GetNodeName() const { return GetNodeNames().front(); }
 
     void AddEventValue(float value) { m_value = value; }
     float GetEventValue() const { return m_value; }
 
     void SetRange(const float a, const float b){
       if (a < b) {
         m_var_range = std::make_pair(a, b);
       } else {
         m_var_range = std::make_pair(b, a);
       }
     }
     void SetMin(float min) { m_var_range.first = min; }
     void SetMax(float max) { m_var_range.second = max; }
     void SetAbsRange(float a) { m_var_range = std::make_pair(-1.0*std::fabs(a), std::fabs(a)); }
     std::pair<float, float> GetRange() const { return m_var_range; }
     float GetMin() const { return m_var_range.first; }
     float GetMax() const { return m_var_range.second; }
 
     std::string GetRangeCutString(const std::string & var="X") const {
       std::string cut_string = Name() + ": ";
       if (std::isnan(GetMin()) && std::isnan(GetMax())) {
         cut_string += "No range set";
       } else if (!std::isnan(GetMin()) && std::isnan(GetMax())) {
         cut_string += var + " > " + std::to_string(GetMin());
       } else if (std::isnan(GetMin()) && !std::isnan(GetMax())) {
         cut_string += var + " < " + std::to_string(GetMax());
       } else if (!std::isnan(GetMin()) && !std::isnan(GetMax())) {
         if (std::fabs(GetMin()) == std::fabs(GetMax())) {
           cut_string += " |" + var + "| < " + std::to_string(std::fabs(GetMin()));
         } else {
           cut_string += std::to_string(GetMin()) + " < " + var + " < " + std::to_string(GetMax());
         }
       }
 
       return cut_string;
     }
 
     bool EvaluateRangeCut(const float x) {
       if (std::isnan(GetMin()) && std::isnan(GetMax())) {
         return true;
       } else if (!std::isnan(GetMin()) && std::isnan(GetMax())) {
         return x > GetMin();
       } else if (std::isnan(GetMin()) && !std::isnan(GetMax())) {
         return x < GetMax();
       } else if (!std::isnan(GetMin()) && !std::isnan(GetMax())) {
         return (x > GetMin() && x < GetMax());
       }
       return false;
     }
 
   protected:
     
     EventCut( const std::string &name = "EventCut" ) :  m_name(name) {}
 
   private:
 
     int m_verbosity{0};
 
     std::string m_name {""};
     bool m_result{false};
 
     std::vector<std::string> m_node_names {};
 
     // only used for cuts that need a range
     std::pair<float, float> m_var_range{NAN, NAN};
 
     float m_value{ NAN };
     std::string m_note{""};
 };
 
 #endif // _EVENTCUT_H_