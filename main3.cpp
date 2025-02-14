#include <cctype>
#include <cassert>
#include <limits>
#include <utility>
#include <algorithm>
#include <iostream>
#include <array>
#include <vector>
#include <unordered_set>
#include <memory>
#include <string>
#include <string_view>
#include <sstream>
#include <charconv>
  
using namespace std::literals::string_view_literals;

enum class LinkType : size_t {
    INVALID_ALTERNATION_RHS_EMPTY,

    NONE,
    CONCATENATION,
    ALTERNATION,
};

enum class OpType : size_t {
    INVALID_REPETITION_UNMATCHED,
    INVALID_REPETITION_CHAR,
    INVALID_REPETITION_COMMA_MORE_THAN_ONE,
    INVALID_REPETITION_N_AND_M_MISSING,
    INVALID_REPETITION_M_LTE_N,

    NONE,
    CONCATENATION,
    ALTERNATION,
    ONE,
    ZERO_OR_ONE,
    ZERO_OR_MORE,
    ONE_OR_MORE,
    N,
    N_OR_MORE,
    M_OR_LESS,
    N_M,
};

enum class GroupType : size_t {
    INVALID_GROUP_UNMATCHED,
    INVALID_GROUP_START,
    INVALID_BRACKET_UNMATCHED,
    INVALID_BRACKET_INNER_UNMATCHED,
    INVALID_BRACKET_INNER_START,
    INVALID_BRACKET_INNER,
    INVALID_BRACKET_RANGE,
    INVALID_BRACKET_COLLATE,
    INVALID_BRACKET_EQUIVALENCE,
    INVALID_BRACKET_HEX,
    INVALID_BRACKET_ESCAPE,

    EMPTY,
    IMPLICIT,
    CAPTURE,
    NONCAPTURE,
    BRACKET,
    BRACKET_NEGATED,
};

struct Expr;

struct Equation {
    std::string text;
    std::vector<Expr*> traversed;
};

struct Expr {
    GroupType group_type;
    OpType op_type;
    LinkType link_type;
    std::string_view group;
    std::string_view op;
    std::string_view link;
    std::vector<std::unique_ptr<Expr>> children;
    std::vector<Equation> eqs;
    size_t n;
    size_t m;

    Expr(GroupType group_type = GroupType::IMPLICIT,
         OpType op_type = OpType::ONE,
         LinkType link_type = LinkType::NONE,
         std::string_view group = "",
         std::string_view op = "",
         std::string_view link = "",
         std::vector<std::unique_ptr<Expr>> children = {},
         size_t n = 0,
         size_t m = 0)
    : group_type(group_type), op_type(op_type), link_type(link_type), group(group), op(op), link(link), children(std::move(children)), n(n), m(m) {}
};

std::string_view scan_number(std::string_view input) {
    size_t max = input.size();
    size_t i = 0;
    for (; i < max; i++) {
        if (!std::isdigit(input[i])) break;
    }
    return input.substr(0, i);
}

OpType check_repetition(std::string_view input) {
    
    // Check that the repetition has an end brace
    auto end_pos = input.find('}');
    if (end_pos == std::string::npos) return OpType::INVALID_REPETITION_UNMATCHED;

    // Check that the repetition content (not including braces) has valid characters
    auto content = input.substr(1, end_pos - 1);
    auto is_valid = [](char c) { return std::isdigit(c) || c == ','; };
    if (!std::all_of(content.begin(), content.end(), is_valid)) return OpType::INVALID_REPETITION_CHAR;

    // Check that there is at most one comma
    if (std::count(content.begin(), content.end(), ',') > 1) return OpType::INVALID_REPETITION_COMMA_MORE_THAN_ONE;

    // If a comma is not found, then the repetition is {n} type
    auto comma_pos = content.find(',');
    if (comma_pos == std::string_view::npos) return OpType::N;

    // if a comma is found at the end, then the repetition is {n,} type
    if (comma_pos == content.size() - 1) return OpType::N_OR_MORE;

    // If a comma is found at the start, then the repetition is {,m} type
    if (comma_pos == 0) return OpType::M_OR_LESS;

    // Otherwise, the repetition is {n,m} type, so check that n < m
    size_t comma_size = 1;
    size_t after_comma_pos = comma_pos + comma_size;
    auto n = scan_number(content);
    auto m = scan_number(content.substr(after_comma_pos));

    size_t nval;
    size_t mval;
    std::from_chars(n.data(), n.data() + n.size(), nval);
    std::from_chars(m.data(), m.data() + m.size(), mval);
    if (nval < mval) return OpType::N_M;
    else return OpType::INVALID_REPETITION_M_LTE_N;
}

OpType check_op(std::string_view input) {
    if (input.empty()) return OpType::ONE;
    auto c = input[0];
    switch (c) {
        case '?': return OpType::ZERO_OR_ONE;
        case '*': return OpType::ZERO_OR_MORE;
        case '+': return OpType::ONE_OR_MORE;
        case '{': return check_repetition(input);
    }
    return OpType::NONE;
}

LinkType check_link(std::string_view input) {
    if (input.empty()) return LinkType::NONE;
    auto c = input[0];
    if (c == '|') {
        if (!(input.size() > 1)) return LinkType::INVALID_ALTERNATION_RHS_EMPTY;
        else return LinkType::ALTERNATION;
    }
    return LinkType::CONCATENATION;
}

bool is_valid_repetition(OpType op_type) {
    return op_type == OpType::N ||
           op_type == OpType::N_OR_MORE ||
           op_type == OpType::M_OR_LESS ||
           op_type == OpType::N_M;
}

bool is_valid_op(OpType op_type) {
    return op_type == OpType::NONE ||
           op_type == OpType::ONE ||
           op_type == OpType::ZERO_OR_ONE ||
           op_type == OpType::ZERO_OR_MORE ||
           op_type == OpType::ONE_OR_MORE ||
           is_valid_repetition(op_type);
}

bool is_valid_link(LinkType link_type) {
    return link_type == LinkType::CONCATENATION ||
           link_type == LinkType::ALTERNATION;
}

void set_range(Expr& expr) {
    if (expr.op.empty()) { expr.n = 0; expr.m = 0; return; }
    auto type = expr.op_type;
    if (!is_valid_op(type)) return;
    auto op = expr.op.substr(1, expr.op.size() - 1);
    auto max = std::numeric_limits<size_t>::max();
    switch (type) {
        case OpType::ZERO_OR_ONE: {
            expr.n = 0;
            expr.m = 1;
            break; 
        }
        case OpType::ZERO_OR_MORE: {
            expr.n = 0;
            expr.m = max;
            break; 
        }
        case OpType::ONE_OR_MORE: {
            expr.n = 1;
            expr.m = max;
            break; 
        }
        case OpType::N: {
            auto n = op;
            std::from_chars(n.data(), n.data() + n.size(), expr.n);
            expr.m = expr.n;
            break;
        }
        case OpType::N_OR_MORE: {
            size_t comma_size = 1;
            auto n = op.substr(0, op.size() - comma_size);
            std::from_chars(n.data(), n.data() + n.size(), expr.n);
            expr.m = max;
            break;
        }
        case OpType::M_OR_LESS: {
            size_t comma_size = 1;
            auto m = op.substr(comma_size);
            expr.n = 0;
            std::from_chars(m.data(), m.data() + m.size(), expr.m);
            break;
        }
        case OpType::N_M: {
            size_t comma_pos = op.find(',');
            auto n = op.substr(0, comma_pos);
            auto m = op.substr(comma_pos + 1);
            std::from_chars(n.data(), n.data() + n.size(), expr.n);
            std::from_chars(m.data(), m.data() + m.size(), expr.m);
            break;
        }
    }
}

std::string_view scan_op(std::string_view input, OpType& op_type) {
    op_type = check_op(input);
    if (op_type == OpType::ZERO_OR_ONE ||
        op_type == OpType::ZERO_OR_MORE ||
        op_type == OpType::ONE_OR_MORE) {
        return input.substr(0, 1);
    }
    else if (is_valid_repetition(op_type)) {
        return input.substr(0, input.find('}') + 1);
    }
    return "";
}

std::string_view scan_link(std::string_view input, LinkType& link_type) {
    link_type = check_link(input);
    if (link_type == LinkType::ALTERNATION) return input.substr(0, 1);
    else return "";
}

bool is_valid_group(GroupType group_type) {
    return group_type == GroupType::IMPLICIT ||
           group_type == GroupType::CAPTURE ||
           group_type == GroupType::NONCAPTURE ||
           group_type == GroupType::BRACKET ||
           group_type == GroupType::BRACKET_NEGATED;
}

GroupType check_group(std::string_view input) {
    if (input.empty()) return GroupType::EMPTY;

    auto c1 = input[0];
    auto c2 = input.size() > 1 ? input[1] : '\0';
    auto c3 = input.size() > 2 ? input[2] : '\0';

    constexpr std::string_view invalid_start = "{|?*+";
    if (invalid_start.contains(c1)) return GroupType::INVALID_GROUP_START;

    if (c1 == '(') {
        return (c2 == '?' && c3 == ':') ? GroupType::NONCAPTURE : GroupType::CAPTURE;
    }
    if (c1 == '[') {
        return (c2 == '^') ? GroupType::BRACKET_NEGATED : GroupType::BRACKET;
    }
    return GroupType::IMPLICIT;
}

// TODO: how to validate a bracketed collation against the current locale?
bool is_valid_bracket_collation(std::string_view input) {
    return true;
}

// TODO: how to validate a bracketed equivalence against the current locale?
bool is_valid_bracket_equivalence(std::string_view input) {
    return true;
}

bool is_valid_bracket_class(std::string_view input) {
    constexpr auto classes = std::array {
        "[:upper:]"sv,
        "[:lower:]"sv,
        "[:alpha:]"sv,
        "[:digit:]"sv,
        "[:xdigit:]"sv,
        "[:alnum:]"sv,
        "[:punct:]"sv,
        "[:blank:]"sv,
        "[:space:]"sv,
        "[:cntrl:]"sv,
        "[:graph:]"sv,
        "[:print:]"sv,
    };
    return std::any_of(classes.begin(), classes.end(), [=](std::string_view sv) {
        return input.starts_with(sv);
    });
}

char char_to_escaped_char(char c) {
    switch (c) {
        case 't': return '\t';
        case 'n': return '\n';
        case 'r': return '\r';
        case 'f': return '\f';
        case 'v': return '\v';
        case 'a': return '\a';
        case 'b': return '\b';
        case '\\': return '\\';
        case '\'': return '\'';
        case '\"': return '\"';
        default: return c;
    }
}

std::string_view scan_bracket_char(std::string_view input, GroupType& type) {
    if (input.empty()) return "";

    char c = input[0];
    bool escaped_char = (c == '\\');
    if (!escaped_char) return input.substr(0, 1);

    size_t escape_size = 2;
    bool escape_fits = (input.size()) >= (escape_size);
    auto escape = (escape_fits) ? (input.substr(0, escape_size)) : ("");
    if (escape_fits && escape != "\\x") return input.substr(0, 2);
    if (!escape_fits) {
        type = GroupType::INVALID_BRACKET_ESCAPE;
        return "";
    }

    size_t skip_size = 2;
    size_t hex_size = 2;
    bool hex_fits = (input.size()) >= (skip_size + hex_size);
    if (hex_fits) {
        auto hex = input.substr(skip_size, hex_size);
        bool hex_valid = std::all_of(hex.begin(), hex.end(), ::isxdigit);
        if (!hex_valid) {
            type = GroupType::INVALID_BRACKET_HEX;
            return "";
        }
        return input.substr(0, skip_size + hex_size);
    }
    else {
        type = GroupType::INVALID_BRACKET_HEX;
        return "";
    }
}

/*
 * precondition: input is well formed (c or \c or \xcc) or ""
 */
char eval_bracket_char(std::string_view input) {
    if (input.empty()) return '\0';
    char c;
    size_t escape_size = 2;
    size_t escaped_hex_size = 4;
    size_t hex_radix = 16;
    if (input.substr(0, escape_size) == "\\x") {
        std::from_chars(input.data() + escape_size, input.data() + escaped_hex_size, c, hex_radix);
        return c;
    }
    if (input[0] == '\\') {
        return char_to_escaped_char(input[1]);
    }
    return input[0];
}

/*
 * precondition: a '[' is at input[0]
 */
std::string_view scan_bracket_inner(std::string_view input, GroupType& type) {
    if (input.size() < 2) {
        type = GroupType::INVALID_BRACKET_INNER_START;
        return "";
    }
    std::string_view valid_second = ":.=";
    auto open = input.substr(0, 2);
    auto second = open[1];
    if (!valid_second.contains(second)) {
        type = GroupType::INVALID_BRACKET_INNER_START;
        return "";
    }
    auto end_it = std::adjacent_find(input.begin(), input.end(), [=](char a, char b) { return (a == second) && (b == ']'); });
    bool closed = (end_it != input.end());
    if (!closed) {
        type = GroupType::INVALID_BRACKET_INNER_UNMATCHED;
        return "";
    }
    else end_it++;
    auto size = std::distance(input.begin(), end_it + 1);
    auto inner = input.substr(0, size);
    if (second == '.' && is_valid_bracket_collation(inner) ||
        second == '=' && is_valid_bracket_equivalence(inner) ||
        second == ':' && is_valid_bracket_class(inner)) {
        return inner;
    }
    else {
       type = GroupType::INVALID_BRACKET_INNER;
       return "";
    }
}

std::string_view scan_bracket(std::string_view input, GroupType& type) {
    constexpr std::string_view special = ":.=";
    auto end_it = std::adjacent_find(input.begin(), input.end(), [=](char a, char b) { return !special.contains(a) && b == ']'; });
    bool closed = (end_it != input.end());
    if (!closed) {
        type = GroupType::INVALID_BRACKET_UNMATCHED;
        return "";
    }
    else end_it++;
    size_t end_pos = std::distance(input.begin(), end_it);
    auto content = input.substr(1 + (input[1] == '^'), end_pos - 1);
    size_t i = 0;
    size_t max = content.size();
    while (i < max) {
        char c = content[i];
        if (c == '[') {
            auto inner = scan_bracket_inner(content, type);
            if (inner.empty()) return "";
            else i += inner.size();
        }
        else {
            auto sub1 = content.substr(i);
            auto sv_c1 = scan_bracket_char(sub1, type);
            char c1 = eval_bracket_char(sv_c1);
            if (c1 == '\0') return "";
            size_t dash_size = 1;
            size_t maybe_dash_pos = (i + sv_c1.size());
            bool dash_fits = (maybe_dash_pos) < (max);
            char maybe_dash = (dash_fits) ? (content[maybe_dash_pos]) : ('\0');
            if (maybe_dash == '-') {
                size_t after_dash_pos = (i + sv_c1.size() + dash_size);
                auto sub2 = content.substr(after_dash_pos);
                auto sv_c2 = scan_bracket_char(sub2, type);
                char c2 = eval_bracket_char(sv_c2);
                if (c2 == '\0') return "";
                if (!(c1 < c2)) {
                    type = GroupType::INVALID_BRACKET_RANGE;
                    return "";
                }
                i += (sv_c1.size() + dash_size + sv_c2.size());
            }
            else i += sv_c1.size();
        }
    }
    return input.substr(0, end_pos + 1);
}

// "ab+" -> "a"
// "ab(" -> "ab"
// "c+" -> "c"
// "c(" -> "c"
std::string_view scan_group_implicit(std::string_view input) {
    constexpr std::string_view end = "(){|?*+";
    constexpr std::string_view split = "?*+{";

    size_t pos_find = input.find_first_of(end);
    if (pos_find == std::string_view::npos) return input;

    auto c = input[pos_find];
    size_t pos_end = (split.contains(c) && pos_find > 1) ? (pos_find - 1) : (pos_find);
    return input.substr(0, pos_end);
}

// "a(..." -> "a"
// "(a)..." -> "(a)"
std::string_view scan_group_explicit(std::string_view input, GroupType& type) {
    int balance = 0;
    size_t i = 0;
    for (; i < input.size(); i++) {
        if (input[i] == '(') {
            balance++;
        }
        else if (input[i] == ')') {
            balance--;
        }
        if (balance == 0) break;
    }
    if (balance != 0) {
        type = GroupType::INVALID_GROUP_UNMATCHED;
        return "";
    }
    return input.substr(0, i + 1);
}

bool is_wrapped(GroupType type) {
    return type == GroupType::CAPTURE ||
           type == GroupType::NONCAPTURE;
}

bool is_bracketed(GroupType type) {
    return type == GroupType::BRACKET ||
           type == GroupType::BRACKET_NEGATED;
}

std::string_view scan_group(std::string_view input, GroupType& type) {
    type = check_group(input);
    bool bracketed = is_bracketed(type);
    if (bracketed) return scan_bracket(input, type);
    bool wrapped = is_wrapped(type);
    if (wrapped) return scan_group_explicit(input, type);
    else return scan_group_implicit(input);
}

// "(a)" -> "a"
std::string_view unwrap_group(std::string_view input) {
    if (input.empty() || input.size() < 2) return "";
    return input.substr(1, input.size() - 2);
}

std::unique_ptr<Expr> parse(std::string_view input) {
    if (input.empty()) return nullptr;

    // Parsing begins here, nests down recursively
    auto root = std::make_unique<Expr>(GroupType::IMPLICIT, OpType::ONE, LinkType::NONE, input, "", "", std::vector<std::unique_ptr<Expr>>{}, 1, 1);

    GroupType group_type;
    auto scan = scan_group(input, group_type);
    bool wrapped = is_wrapped(group_type);

    // Handle final nesting
    if ((wrapped ? unwrap_group(scan) : scan).size() == input.size()) {
        root->group_type = group_type;
        return std::move(root);
    }

    // TODO: Handle group parse errors
    if (!is_valid_group(group_type)) {}

    auto rest = input.substr(scan.size());

    // Handle increasing nesting and left-to-right group traversal
    while (scan.size() > 0) {

        // Parse the left hand group
        auto maybe_unwrap = !wrapped ? scan : unwrap_group(scan);
        auto lhs = parse(maybe_unwrap);

        // Scan the operation applied to lhs
        lhs->group = scan;
        lhs->group_type = group_type;
        lhs->op = scan_op(rest, lhs->op_type);
        lhs->link = scan_link(rest.substr(lhs->op.size()), lhs->link_type);
        set_range(*lhs);

        // TODO: Handle op and link parse errors
        if (!is_valid_op(lhs->op_type)) {}
        if (!is_valid_link(lhs->link_type)) {}

        // Consume the operation and link (concatenation removes nothing)
        rest.remove_prefix(lhs->op.size() + lhs->link.size());

        // Add the group as a child
        root->children.push_back(std::move(lhs));

        // Set up the next loop iteration
        scan = scan_group(rest, group_type);
        wrapped = is_wrapped(group_type);
        rest = rest.substr(scan.size());

        // TODO: Handle group parse errors
        if (!is_valid_group(group_type)) {}
    }

    return std::move(root);
}

bool is_number(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

void merge_eqs(std::vector<Equation>& lhs, const std::vector<Equation>& rhs) {
    for (const auto& eq : rhs) {
        lhs.push_back(eq);
    }
}

std::string_view scan_const(std::string_view input) {
    auto maybe_const = scan_number(input);
    auto rest = input.substr(maybe_const.size());
    return (rest.empty() || rest[0] == '+') ?
           (maybe_const) :
           ("");
}

void cart_prod_pair_sum_eqs(std::vector<Equation>& lhs, const std::vector<Equation>& rhs) {
    std::vector<Equation> new_lhs;
    for (auto& left : lhs) {
        auto left_const = scan_const(left.text);
        auto left_expr = left.text.substr(left_const.size());
        size_t left_val = 0;
        std::from_chars(left_const.data(), left_const.data() + left_const.size(), left_val);
        for (auto& right : rhs) {
            auto right_const = scan_const(right.text);
            auto right_expr = right.text.substr(right_const.size());
            size_t right_val = 0;
            std::from_chars(right_const.data(), right_const.data() + right_const.size(), right_val);
            int combined_const = left_val + right_val;
            std::string final_expr = ((combined_const != 0) ?
                                         (std::to_string(combined_const)) :
                                         ("")) +
                                     ((!left_expr.empty() && (left_expr[0] != '+') && combined_const != 0) ?
                                         ("+") :
                                         ("")) +
                                     std::string(left_expr) +
                                     ((!right_expr.empty() && (right.text[0] != '+')) ?
                                         ("+") :
                                         ("")) +
                                     std::string(right_expr);
            std::vector<Expr*> traversed = left.traversed;
            for (auto* t : right.traversed) {
                traversed.push_back(t);
            }
            new_lhs.emplace_back(std::move(final_expr), std::move(traversed));
        }
    }
    lhs = new_lhs;
}

void scalar_mult_eqs(std::vector<Equation>& eqs, std::string& c) {
    for (auto& eq : eqs) {
        size_t last_pos = 0;
        size_t pos = 0;
        std::string_view term;
        while ((pos = eq.text.find('+', pos)) != std::string::npos) {
            term = eq.text.substr(last_pos, pos);
            if (term[0] == '{') {
                eq.text.insert(last_pos, "1");
                pos++;
            }
            eq.text.insert(pos, c);
            pos += c.length() + 1;
            last_pos = pos;
        }
        term = eq.text.substr(last_pos);
        if (term[0] == '{') {
            eq.text.insert(last_pos, "1");
            pos++;
        }
        eq.text += c;
    }
}

std::string n_m_to_var(size_t n, size_t m, size_t var_count) {
    auto nt = std::to_string(n);
    auto mt = std::to_string(m);
    auto count = std::to_string(var_count);
    return "{x" + count + ":" + nt + "," + mt + "}";
}

void gen_eqs(Expr& expr, size_t& var_count) {
    Expr* first = nullptr;
    LinkType link_type = LinkType::NONE;
    if (expr.children.empty()) {
        auto number = std::to_string(expr.group.size());
        expr.eqs.emplace_back(number, std::vector<Expr*>{&expr});
    }
    else {
        first = expr.children[0].get();
        gen_eqs(*first, var_count);
        merge_eqs(expr.eqs, first->eqs);
        link_type = first->link_type;
    }
    size_t max = expr.children.size();
    for (size_t i = 1; i < max; i++) {
        Expr& ch = *expr.children[i];
        gen_eqs(ch, var_count);
        switch (link_type) {
            case LinkType::ALTERNATION: {
                merge_eqs(expr.eqs, ch.eqs); break;
            }
            case LinkType::CONCATENATION: {
                cart_prod_pair_sum_eqs(expr.eqs, ch.eqs); break;
            }
            case LinkType::NONE: {
                cart_prod_pair_sum_eqs(expr.eqs, ch.eqs); break;
            }
        }
        link_type = ch.link_type;
    }
    if (expr.op_type != OpType::NONE && expr.op_type != OpType::ONE) {
        auto var = n_m_to_var(expr.n, expr.m, var_count++);
        scalar_mult_eqs(expr.eqs, var);
    }
}

// ------------------------------------------------------------------
// Helper: convert generated equation string to solver format.
// Our internal format looks like:
//    "1*{x0:3,3}*{x1:1,1}+2*{x1:0,3}*{x2:2,5}=9"
// We want to convert that to:
//    "1*x0*x1+2*x1*x2=9"
// That is, remove each brace-enclosed bound, leaving only the variable name.
// ------------------------------------------------------------------
std::string format_eq(const Equation& eq, std::string& input) {
    std::string result;
    result.reserve(eq.text.size());
    // We scan the string character by character.
    for (size_t i = 0; i < eq.text.size(); ++i) {
        char c = eq.text[i];
        if (c == '{') {
            // We assume the format is "{<var>:<n>,<m>}"
            // Find the colon and then the closing brace.
            size_t colonPos = eq.text.find(':', i);
            size_t closePos = eq.text.find('}', i);
            if (colonPos == std::string::npos || closePos == std::string::npos) {
                // If something is wrong, copy the rest and break.
                result.append(eq.text.substr(i));
                break;
            }
            // Append the substring from after '{' up to the colon.
            result.append("*");
            result.append(eq.text.substr(i + 1, colonPos - i - 1));
            // Skip ahead past the closing brace.
            i = closePos;
        } else {
            result.push_back(c);
        }
    }
    result.append("=" + std::to_string(input.size()));
    return result;
}

// ------------------------------------------------------------------
// Helper: call the Python solver program (solver.py)
// with a given equation string (in solver format) as an argument,
// capture its output and return it as a string.
// ------------------------------------------------------------------
std::string solve_eq(const std::string &equation) {
    // Construct command. We enclose the equation in quotes.
    std::string command = "python solver.py \"" + equation + "\"";
    FILE *pipe = popen(command.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Unable to run solver.py." << std::endl;
        return "";
    }
    const int bufSize = 256;
    char buffer[bufSize];
    std::stringstream output;
    while (fgets(buffer, bufSize, pipe) != nullptr) {
        output << buffer;
    }
    pclose(pipe);
    return output.str();
}

// Function to print parsed expressions for debugging
void print_expr(const Expr& expr, int indent = 0) {
    std::string padding(indent * 2, ' ');
    std::cout << padding
              << "Expr:" << std::endl
              << ", Group = \"" << expr.group
              << "\", Op = \"" << expr.op
              << "\", Link = \"" << expr.link
              << "\", Op String = \"" << expr.op << "\""
              << ", n = " << expr.n
              << ", m = " << expr.m
              << ", m = " << expr.m;

    for (auto& eq : expr.eqs) {
        std::cout << padding
                  << eq.text << std::endl;
    }

    for (const auto& child : expr.children) {
        print_expr(*child, indent + 1);
    }
}

void test_parse() {
    std::cout << "Running parse tests...\n";

    /*
    // Test 1: Simple character
    {
        auto expr = parse("a");
        assert(expr && expr->group == "a");
        assert(expr->op_type == OpType::ONE);
    }

    // Test 2: Simple repetition {3}
    {
        auto expr = parse("a{3}");
        assert(expr && expr->children.size() == 1);
        assert(expr->children[0]->op_type == OpType::N);
        assert(expr->children[0]->n == 3);
        assert(expr->children[0]->m == 3);
    }

    // Test 3: Repetition with range {2,5}
    {
        auto expr = parse("a{2,5}");
        assert(expr && expr->children.size() == 1);
        assert(expr->children[0]->op_type == OpType::N_M);
        assert(expr->children[0]->n == 2);
        assert(expr->children[0]->m == 5);
    }

    // Test 4: Group
    {
        auto expr = parse("(ab)");
        assert(expr && expr->group == "(ab)");
        assert(expr->group_type == GroupType::CAPTURE);
    }

    // Test 5: Alternation
    {
        auto expr = parse("a|b");
        assert(expr && expr->group == "a|b" && expr->op_type == OpType::ONE && expr->children.size() == 2);
        assert(expr->children[0]->group == "a");
        assert(expr->children[1]->group == "b");
        assert(expr->children[0]->op_type == OpType::ALTERNATION);
        assert(expr->children[1]->op_type == OpType::ONE);
    }

    // Test 6: Group with operator 
    {
        auto expr = parse("(ab)*");
        assert(expr && expr->group == "(ab)*" && expr->op_type == OpType::ONE && expr->children.size() == 1);
        assert(expr->children[0]->group == "(ab)");
        assert(expr->children[0]->op_type == OpType::ZERO_OR_MORE);
    }
    */

    // Test 7: Nested groups
    {
        auto expr = parse("((ab)c)+");
        assert(expr && expr->group == "((ab)c)+" && expr->op_type == OpType::ONE && expr->children.size() == 1);
    }

    // Test 8: Some bracket stuff
    {
        auto expr = parse("[[:alnum:]\\x61-\\x7A\t]");
    }

}

int main() {
    /*
    test_parse();
    std::string input = "((a|bc|d{1,5})(e|fg|h{2,3})){4,6}";
    std::cout << "parsing " << input << std::endl;
    auto expr = parse(input);
    print_expr(*expr, 0);
    size_t count = 0;
    auto eqs = gen_eqs(*expr, count);
    std::cout << "equations" << std::endl;
    for (auto& eq : eqs) std::cout << eq << std::endl;
    */
    std::string regex;
    std::cout << "Enter a regex: " << std::endl;
    std::getline(std::cin, regex);

    std::string input;
    std::cout << "Enter an input: " << std::endl;
    std::getline(std::cin, input);
    
    // Step 1. Parse the input using your parse() function.
    // (Assume that parse returns a std::unique_ptr<Expr>.)
    std::unique_ptr<Expr> expr = parse(regex);
    if (!expr) {
        std::cerr << "Error: Parsing failed." << std::endl;
        return 1;
    }
    
    // Step 2. Call gen_eqs() with a counter variable.
    size_t var_count = 0;
    // gen_eqs() returns a vector of strings representing equations in your internal format.
    gen_eqs(*expr, var_count);
    if (expr->eqs.empty()) {
        std::cerr << "Error: gen_eqs produced no equations." << std::endl;
        return 1;
    }

    print_expr(*expr);
    
    // Debug: print the generated equations.
    std::cout << "Generated equations (internal format):" << std::endl;
    for (const auto &eq : expr->eqs) {
        std::cout << eq.text << std::endl;
    }
    
    // Step 3. Convert each generated equation into the format expected by the Python solver.
    std::vector<std::string> formatted_eqs;
    for (const auto &eq : expr->eqs) {
        std::string feq = format_eq(eq, input);
        formatted_eqs.push_back(feq);
    }
    
    // Debug: print the formatted equations.
    std::cout << "\nFormatted equations (solver format):" << std::endl;
    for (const auto &eq : formatted_eqs) {
        std::cout << eq << std::endl;
    }
    
    // Step 4. For each converted equation, call the Python solver and output its result.
    std::cout << "\nSolutions from Python solver:" << std::endl;
    for (const auto &eq : formatted_eqs) {
        std::cout << "For equation: " << eq << std::endl;
        std::string solverOutput = solve_eq(eq);
        std::cout << solverOutput << std::endl;
    }
    return 0;
}

