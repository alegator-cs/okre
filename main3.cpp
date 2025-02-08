#include <cctype>
#include <cassert>
#include <limits>
#include <utility>
#include <algorithm>
#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include <string>
#include <string_view>
#include <charconv>
  
using namespace std::literals::string_view_literals;

enum class OpType : size_t {
    UNREACHABLE,
    INVALID_REPETITION_UNMATCHED,
    INVALID_REPETITION_CHAR,
    INVALID_REPETITION_COMMA_MORE_THAN_ONE,
    INVALID_REPETITION_N_AND_M_MISSING,
    INVALID_REPETITION_M_LTE_N,

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

struct Expr {
    GroupType group_type;
    OpType op_type;
    std::string_view group;
    std::string_view op;
    std::vector<std::unique_ptr<Expr>> children;
    size_t n;
    size_t m;

    Expr(GroupType group_type = GroupType::IMPLICIT,
         OpType op_type = OpType::ONE,
         std::string_view group = "",
         std::string_view op = "",
         std::vector<std::unique_ptr<Expr>> children = {},
         size_t n = 0,
         size_t m = 0)
    : group_type(group_type), op_type(op_type), group(group), op(op), children(std::move(children)), n(n), m(m) {}
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
        case '|': return OpType::ALTERNATION;
        case '{': return check_repetition(input);
    }
    return OpType::CONCATENATION;
}

bool is_valid_repetition(OpType op_type) {
    return op_type == OpType::N ||
           op_type == OpType::N_OR_MORE ||
           op_type == OpType::M_OR_LESS ||
           op_type == OpType::N_M;
}

bool is_valid_op(OpType op_type) {
    return op_type == OpType::ONE ||
           op_type == OpType::ZERO_OR_ONE ||
           op_type == OpType::ZERO_OR_MORE ||
           op_type == OpType::ONE_OR_MORE ||
           op_type == OpType::ALTERNATION ||
           op_type == OpType::CONCATENATION ||
           is_valid_repetition(op_type);
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
        case OpType::ALTERNATION: {
            expr.n = 0;
            expr.m = 1;
            break; 
        }
        case OpType::CONCATENATION: {
            expr.n = 1;
            expr.m = 1;
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
        op_type == OpType::ONE_OR_MORE ||
        op_type == OpType::ALTERNATION) {
        return input.substr(0, 1);
    }
    else if (op_type != OpType::CONCATENATION) {
        if (is_valid_repetition(op_type)) {
            return input.substr(0, input.find('}') + 1);
        }
    }
    return "";
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
    size_t char_size = 1;
    size_t escape_size = 2;
    size_t hex_size = 2;
    bool escape_fits = (input.size()) >= (escape_size);
    if (escape_fits) {
        auto maybe_escape_hex = input.substr(0, escape_size);
        if (maybe_escape_hex == "\\x") {
            size_t min_size_for_hex = 4;
            bool hex_fits = (input.size()) >= (min_size_for_hex);
            auto maybe_hex = input.substr(escape_size, hex_size);
            bool hex_valid = (hex_fits) && (std::isxdigit(maybe_hex[0])) && (std::isxdigit(maybe_hex[1]));
            if (!hex_valid) {
                type = GroupType::INVALID_BRACKET_HEX;
                return "";
            }
            return input.substr(0, escape_size + hex_size);
        }
        else if (c == '\\') {
            return input.substr(0, escape_size);
        }
    }
    else if (c == '\\') {
        type = GroupType::INVALID_BRACKET_ESCAPE;
        return "";
    }
    return input.substr(0, char_size);
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
    auto unwrapped = input.substr(1 + (input[1] == '^'), end_pos - 1);
    size_t i = 0;
    size_t max = unwrapped.size();
    while (i < max) {
        char c = unwrapped[i];
        if (c == '[') {
            auto inner = scan_bracket_inner(unwrapped, type);
            if (inner.empty()) return "";
            else i += inner.size();
        }
        else {
            auto sub1 = unwrapped.substr(i);
            auto sv_c1 = scan_bracket_char(sub1, type);
            char c1 = eval_bracket_char(sv_c1);
            if (c1 == '\0') return "";
            size_t dash_size = 1;
            size_t maybe_dash_pos = (i + sv_c1.size());
            bool dash_fits = (maybe_dash_pos) < (max);
            char maybe_dash = (dash_fits) ? (unwrapped[maybe_dash_pos]) : ('\0');
            if (maybe_dash == '-') {
                size_t after_dash_pos = (i + sv_c1.size() + dash_size);
                auto sub2 = unwrapped.substr(after_dash_pos);
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

// "abc+" -> "ab"
// "abc(" -> "abc"
// "c+" -> ""
// "c(" -> "c"
std::string_view scan_group_implicit(std::string_view input) {
    constexpr std::string_view end = "(){|?*+";
    constexpr std::string_view split = "?*+";

    size_t pos_find = input.find_first_of(end);
    if (pos_find == std::string_view::npos) return input;

    auto c = input[pos_find];
    size_t pos_end = (split.contains(c) && pos_find > 0) ? (pos_find - 1) : (pos_find);
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

std::string_view scan_group(std::string_view input, GroupType& type) {
    type = check_group(input);
    if (type == GroupType::BRACKET || type == GroupType::BRACKET_NEGATED) return scan_bracket(input, type);
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
    auto root = std::make_unique<Expr>(GroupType::IMPLICIT, OpType::ONE, input, "", std::vector<std::unique_ptr<Expr>>{}, 1, 1);

    GroupType group_type;
    auto scan = scan_group(input, group_type);

    // Handle final nesting
    if (scan.size() == input.size()) {
        root->group_type = group_type;
        return std::move(root);
    }

    // TODO: Handle group parse errors
    if (!is_valid_group(group_type)) {}

    bool wrapped = is_wrapped(group_type);
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
        set_range(*lhs);

        // TODO: Handle op parse errors
        if (!is_valid_op(lhs->op_type)) {}

        // Consume the operation (concatenation removes nothing)
        rest.remove_prefix(lhs->op.size());

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

// Function to print parsed expressions for debugging
void print_expr(const Expr& expr, int indent = 0) {
    std::string padding(indent * 2, ' ');
    std::cout << padding << "Expr: Group = " << static_cast<size_t>(expr.group_type)
              << ", Op = " << static_cast<size_t>(expr.op_type)
              << ", n = " << expr.n
              << ", m = " << expr.m
              << ", Group String = \"" << expr.group
              << "\", Op String = \"" << expr.op << "\"\n";

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

    std::cout << "All tests passed!\n";
}

int main() {
    test_parse();
    return 0;
}

