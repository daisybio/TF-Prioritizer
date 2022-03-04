const textDownloadPrefix = "data:text/plain;charset=utf-8,";

function toggleAccordion(id) {
    let panel = document.getElementById(id);
    if (panel.style.display === "grid") {
        panel.style.display = "none";

        let info = document.getElementById(id + "-info");
        if (info != null) {
            info.style.display = "none";
        }
    } else {
        panel.style.display = "grid";
    }
}

function openImageInTab(id) {
    let image = document.getElementById(id);
    window.open(image.src);
}

function openModal(modal_id) {
    let modal = document.getElementById(modal_id);
    modal.style.display = "block";
}

function closeModal(modal_id) {
    let modal = document.getElementById(modal_id);
    modal.style.display = "none";
}

function tableMouseOver(gene, col, row) {
    colHead = document.getElementById(gene + "-col-" + col);
    rowHead = document.getElementById(gene + "-row-" + row);

    if (colHead !== null) {
        colHead.style.backgroundColor = "rgba(var(--black), 0.1)";
    }
    if (rowHead !== null) {
        rowHead.style.backgroundColor = "rgba(var(--black), 0.1)";
    }
}

function tableMouseOut(gene, col, row) {
    colHead = document.getElementById(gene + "-col-" + col);
    rowHead = document.getElementById(gene + "-row-" + row);

    if (colHead !== null) {
        colHead.style.backgroundColor = "initial";
    }
    if (rowHead !== null) {
        rowHead.style.backgroundColor = "initial";
    }
}

function getAllCombinationEntries(combinations) {
    if (Array.isArray(combinations)) {
        return new Set(combinations);
    } else {
        let childrenSet = new Set();
        for (let i = 0; i < Object.keys(combinations).length; i++) {
            let addSet = getAllCombinationEntries(combinations[Object.keys(combinations)[i]]);
            addSet.forEach(entry => childrenSet.add(entry));
        }
        return childrenSet;
    }
}

function init_filterOptions(id) {
    let combinations = window[id + "CombinationsUnfiltered"];
    let filterOptions = document.getElementsByClassName("filterOption " + id);

    if (filterOptions.length > 0) {
        let foundFilterOptions = new Set();
        let regexTemplate = "^[0-9]+_{PREFIX}.+"
        let entries = getAllCombinationEntries(combinations);

        let symbolFilterOption;
        for (let i = 0; i < filterOptions.length; i++) {
            if (filterOptions[i].id.endsWith("Symbols")) {
                symbolFilterOption = filterOptions[i];
            }
        }

        entries.forEach(entry => {
            let assigned = false;
            for (let i = 0; i < filterOptions.length; i++) {
                if (filterOptions[i] !== symbolFilterOption) {
                    let regexString = regexTemplate.replace("{PREFIX}", filterOptions[i].value);
                    if (new RegExp(regexTemplate.replace("{PREFIX}", filterOptions[i].value)).test(entry)) {
                        foundFilterOptions.add(filterOptions[i]);
                        assigned = true;
                    }
                }
            }

            if (!assigned) {
                foundFilterOptions.add(symbolFilterOption);
            }
        });

        for (let filterOption of filterOptions) {
            if (foundFilterOptions.has(filterOption)) {
                filterOption.classList.add("active");
                filterOption.addEventListener("click", () => {
                    toggleFilterActive(filterOption, id);
                    window[id + "Combinations"] = getFilteredCombinations(id);
                    init_selection(id, window[id + "Combinations"]);
                });
            } else {
                filterOption.disabled = true;
                filterOption.title = "This gene type is not present in this analysis"
            }
        }

        if (countActiveElements(filterOptions) === 1) {
            getActiveElement(filterOptions).click();
            //filterOptions[0].parentElement.style.display = "none";
        }
    }
}

function init_selection(id) {
    let combinations = window[id + "Combinations"];
    let firstValue;
    let availableValues;

    if (Array.isArray(combinations)) {
        availableValues = combinations;
    } else {
        availableValues = Object.keys(combinations);
    }

    let activeChild = getActiveElement(document.querySelectorAll('[id^="' + id + '-0"]'));
    if (activeChild !== null && availableValues.includes(activeChild.value)) {
        activeChild.click();
    } else {
        let term = id + "-0-" + availableValues[0];
        let firstChild = document.getElementById(term);
        if (firstChild !== null) {
            firstChild.click();
        }
    }
}

function update_selection(source_element, id) {
    let combinations = window[id + "Combinations"];
    let level = parseInt(source_element.id.split("-")[1]);
    let depth = get_depth(combinations);

    let parentElements = [];

    { // Find activated parent nodes
        for (let i = 0; i < level; i++) {
            let nodes = document.querySelectorAll('[id^="' + id + '-' + i + '"]');

            for (let j = 0; j < nodes.length; j++) {
                if (nodes[j].classList.contains("active")) {
                    parentElements.push(nodes[j].value);
                    break;
                }
            }
        }
    }

    let availableCombinations = combinations;

    for (let i = 0; i < parentElements.length; i++) {
        availableCombinations = availableCombinations[parentElements[i]];
    }

    { // Activate clicked element and deactivate siblings
        let siblings = document.querySelectorAll('[id^="' + id + '-' + level + '"]');
        let allowedSiblings;

        if (Array.isArray(availableCombinations)) {
            allowedSiblings = availableCombinations;
        } else {
            allowedSiblings = Object.keys(availableCombinations);
        }

        for (let i = 0; i < siblings.length; i++) {
            siblings[i].classList.remove("active");
            siblings[i].disabled = !allowedSiblings.includes(siblings[i].value);
        }
        source_element.classList.add("active");
    }

    availableCombinations = availableCombinations[source_element.value];

    { // Fill dropdown if necessary
        if (level === depth - 1 && document.getElementById(id + "-dropdown") != null) {
            let dropdownContent = document.getElementById(id + "-dropdown-content");
            dropdownContent.innerHTML = "";

            let searchbox = document.createElement("input");
            searchbox.type = "text";
            searchbox.placeholder = "Search...";
            searchbox.id = id + "-dropdown-search";
            searchbox.style.width = "100%";
            searchbox.addEventListener("keyup", function () {
                filter_dropdown(id)
            });
            searchbox.classList.add("selector")

            dropdownContent.appendChild(searchbox);

            let filterOptions = document.getElementsByClassName("filterOption " + id);
            let regexTemplate = "^[0-9]+_{PREFIX}.+";

            let symbolFilterOption;
            for (let i = 0; i < filterOptions.length; i++) {
                if (filterOptions[i].id.endsWith("Symbols")) {
                    symbolFilterOption = filterOptions[i];
                    break;
                }
            }

            for (let i = 0; i < availableCombinations.length; i++) {
                let option = document.createElement("button");
                let value = availableCombinations[i];

                let keepEntry = true;

                if (filterOptions.length > 0) {
                    let assigned = false;


                    for (let j = 0; j < filterOptions.length; j++) {
                        if (filterOptions[j] !== symbolFilterOption) {
                            let regexString = regexTemplate.replace("{PREFIX}", filterOptions[j].value);
                            if (new RegExp(regexString).test(value)) {
                                keepEntry = filterOptions[j].classList.contains("active");
                                assigned = true;
                                break;
                            }
                        }
                    }
                    if (!assigned) {
                        keepEntry = symbolFilterOption.classList.contains("active");
                    }
                }

                if (keepEntry) {
                    option.value = value;
                    option.id = id + "-" + 2 + "-" + value;
                    option.classList.add("selector");

                    let text = value.replace(/\.[^/.]+$/, "");
                    if (/^[0-9]+_*/.test(text)) {
                        text = text.substring(text.split("_")[0].length + 1);
                    }

                    option.textContent = (i + 1) + ". " + text;
                    option.addEventListener("click", function () {
                        update_selection(option, id, combinations);
                    });

                    dropdownContent.appendChild(option);
                }
            }

            let dropdownButton = document.getElementById(id + "-dropdown");
            let dropdownLeftArrow = document.getElementById(id + "-next-option");
            let dropdownRightArrow = document.getElementById(id + "-previous-option");
            let modalLeftArrow = document.getElementById(id + "-modal-leftarrow");
            let modalRightArrow = document.getElementById(id + "-modal-rightarrow");

            dropdownButton.disabled = dropdownLeftArrow.disabled = dropdownRightArrow.disabled = modalLeftArrow.disabled = modalRightArrow.disabled = availableCombinations.length <= 1;
        }
    }

    { // Click first possible child
        if (level < depth) {
            let availableValues = [];

            if (Array.isArray(availableCombinations)) {
                availableValues = availableCombinations;
            } else {
                availableValues = Array.from(Object.keys(availableCombinations));
            }

            let currentlyActiveChild = getActiveElement(document.querySelectorAll('[id^="' + id + '-' + (level + 1) + '"]'));

            if (currentlyActiveChild !== null && availableValues.includes(currentlyActiveChild.value)) {
                currentlyActiveChild.click()
            } else {
                let child;
                let found = false;

                for (let i = 0; i < availableValues.length; i++) {
                    let term = id + "-" + (level + 1) + "-" + availableValues[i];
                    child = document.getElementById(term);
                    if (child !== null) {
                        found = true;
                        break;
                    }
                }

                if (found) {
                    child.disabled = false;
                    child.click();
                }
            }
        }
    }

    { // Update image if necessary
        if (level === depth) {
            let image = document.getElementById(id + "-image");
            let modal_image = document.getElementById(id + "-modal-image");

            let source = id + "/";

            for (let i = 0; i < parentElements.length; i++) {
                source += parentElements[i] + "/";
            }
            source += source_element.value;

            modal_image.src = image.src = source;

            let modal_caption = document.getElementById(id + "-modal-caption");
            modal_caption.textContent = source_element.textContent;

            let dropdown = document.getElementById(id + "-dropdown");
            if (dropdown != null) {
                dropdown.textContent = source_element.textContent;
                dropdown.value = source_element.value;
            }
        }
    }
}

function move_lowest_level(id, delta) {
    let combinations = window[id + "Combinations"];
    let options = document.querySelectorAll('[id^="' + id + '-' + get_depth(combinations) + '"]');
    let i;

    for (i = 0; i < options.length; i++) {
        if (options[i].classList.contains("active")) {
            break;
        }
    }
    i += delta;
    if (i >= options.length) {
        i = 0;
    }
    if (i < 0) {
        i = options.length - 1;
    }

    options[i].click();

}

function openPanel(id) {
    let panel = document.getElementById(id);
    panel.style.display = "grid";
}

function toggleDropdown(selection) {
    document.getElementById(selection + "-dropdown-content").classList.toggle("show");
}

function filter_dropdown(selection) {
    let searchbox = document.getElementById(selection + "-dropdown-search");
    let term = searchbox.value.toUpperCase();
    let options = document.getElementById(selection + "-dropdown-content").childNodes;

    for (let i = 0; i < options.length; i++) {
        if (options[i] === searchbox) {
            continue;
        }
        if (options[i].textContent.toUpperCase().indexOf(term) > -1) {
            options[i].style.display = "";
        } else {
            options[i].style.display = "none";
        }
    }
}

function add_dropdown_closing(id) {
    document.addEventListener("click", function (e) {
        let content = document.getElementById(id + "-dropdown-content");
        if (!content.contains(e.target) && e.target !== document.getElementById(id + "-dropdown")) {
            content.classList.remove("show");
        }
    });
}

function get_depth(combinations) {
    let i = 0;

    while (!Array.isArray(combinations)) {
        combinations = combinations[Object.keys(combinations)[0]];
        i++;
    }
    return i;
}

function downloadVisibleTfs() {
    let button = document.getElementById("tfDownload");

    let tfs = document.getElementsByClassName("tf");

    let tfString = "";

    for (let i = 0; i < tfs.length; i++) {
        if (tfs[i].style.display != "none") {
            tfString += tfs[i].id + "\n";
        }
    }

    let fileName = "tfs_";

    if (document.getElementById("searchTFNames").value != "") {
        fileName += "nameFilter_" + document.getElementById("searchTFNames").value;
    } else if (document.getElementById("searchTargetGenes").value != "") {
        fileName += "targetGeneFilter_" + document.getElementById("searchTargetGenes").value;
    } else {
        fileName += "all";
    }

    download(tfString, fileName + ".csv");
}

function filterTfsByName(term) {
    let tfs = document.getElementsByClassName("tf");

    for (let i = 0; i < tfs.length; i++) {
        if (tfs[i].id.toUpperCase().indexOf(term.toUpperCase()) > -1) {
            tfs[i].style.display = "";
        } else {
            tfs[i].style.display = "none";
        }
    }
}

function filterTfsByTargetGene(term) {
    let tfs = document.getElementsByClassName("tf");

    for (let i = 0; i < tfs.length; i++) {
        let varName = "targetGenes" + tfs[i].id;
        let targetGenes = window[varName];
        let found = false;

        for (let j = 0; j < targetGenes.length; j++) {
            if (targetGenes[j].toUpperCase().indexOf(term.toUpperCase()) > -1) {
                found = true;
                break;
            }
        }

        if (found) {
            tfs[i].style.display = "";
        } else {
            tfs[i].style.display = "none";
        }
    }
}

function toggleFilterActive(element, id) {
    let filterOptions = document.getElementsByClassName("filterOption " + id);

    if (element.classList.contains("active")) {
        if (countActiveElements(filterOptions) > 1) {
            element.classList.remove("active");
        }
        if (countActiveElements(filterOptions) === 1) {
            let lastActive = getActiveElement(filterOptions);
            lastActive.style.cursor = "not-allowed";
            lastActive.title = "At least one gene type has to be active";
        }
    } else {
        element.classList.add("active");

        for (let filterOption of filterOptions) {
            if (!filterOption.disabled) {
                filterOption.style.cursor = "pointer";
                filterOption.removeAttribute("title");
            }
        }
    }
}

function getActiveElement(elements) {
    for (let element of elements) {
        if (element.classList.contains("active")) {
            return element;
        }
    }
    return null;
}

function getFilteredCombinations(id) {
    let allCombinations = window[id + "CombinationsUnfiltered"];
    let filterOptions = document.getElementsByClassName("filterOption " + id);
    let regexTemplate = "^[0-9]+_{PREFIX}.+";

    let symbolFilterOption;
    for (let filterOption of filterOptions) {
        if (filterOption.id.endsWith("Symbols")) {
            symbolFilterOption = filterOption;
            break;
        }
    }

    function filterCombination(allCombinations) {
        if (Array.isArray(allCombinations)) {
            let filteredCombinations = [];

            for (let entry of allCombinations) {
                let keep = false;
                for (let filterOption of filterOptions) {
                    if (filterOption !== symbolFilterOption && new RegExp(regexTemplate.replace("{PREFIX}", filterOption.value)).test(entry)) {
                        keep = filterOption.classList.contains("active");
                    }
                }
                if (!keep) {
                    keep = symbolFilterOption.classList.contains("active");
                }
                if (keep) {
                    filteredCombinations.push(entry);
                }
            }

            return filteredCombinations;
        } else {
            let filteredCombinations = {};

            for (let [key, value] of Object.entries(allCombinations)) {
                let values = filterCombination(value);

                if (Array.isArray(values) && values.length > 0) {
                    filteredCombinations[key] = values;
                } else if (Object.keys(values).length > 0) {
                    filteredCombinations[key] = values;
                }
            }

            return filteredCombinations;
        }
    }

    return filterCombination(allCombinations);
}

function download(text, filename) {
    let element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + text);
    element.setAttribute('download', filename);

    element.style.display = 'none';
    document.body.appendChild(element);

    element.click();

    document.body.removeChild(element);
}

function countActiveElements(elements) {
    let activeCount = 0;

    for (let element of elements) {
        if (element.classList.contains("active")) {
            activeCount++;
        }
    }

    return activeCount;
}

function downloadActive(id) {
    let dataCombinations = window[id + "DataCombinations"];
    let level = 0;
    let steps = "";

    while (typeof dataCombinations !== "string") {
        let active = getActiveElement(document.querySelectorAll('[id^=\"' + id + '-' + level + '\"]'));
        steps += active.value + "-";
        dataCombinations = dataCombinations[active.value];
    }

    steps = steps.substring(0, steps.length - 1);

    download(dataCombinations, id + "_" + steps + ".csv");
}