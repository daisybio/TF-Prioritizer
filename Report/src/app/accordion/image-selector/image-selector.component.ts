import {Component, ElementRef, HostListener, Input, OnInit} from '@angular/core';
import {dataPrefix} from "../../../assets/constants";

@Component({
  selector: 'image-selector',
  templateUrl: './image-selector.component.html',
  styleUrls: ['./image-selector.component.css']
})
export class ImageSelectorComponent implements OnInit {
  @Input()
    // @ts-ignore
  data: {};

  @Input()
  visible: boolean = false;

  @Input()
  title: string = "undefined";

  @Input()
  sub: boolean = false;

  @Input()
  information: string | undefined;

  availableOptions: string[][] = [];
  activeOptions: string[] = [];
  allowedOptions: string[][] = [];

  hasDropDown: boolean = false;
  hasData: boolean = false;
  hasPlots: boolean = false;

  isEnabled: boolean = true;

  dropdownVisible: boolean = false;

  imageSource: string = "";

  constructor(private eRef: ElementRef) {
  }

  isAllowed(level: number, option: string) {
    return this.allowedOptions[level].includes(option);
  }

  ngOnInit(): void {
    this.collectLevels(this.data, 0);

    this.hasData = this.availableOptions[this.availableOptions.length - 1].includes("data");
    this.hasPlots = this.availableOptions[this.availableOptions.length - 1].includes("plot");

    this.isEnabled = this.hasData || this.hasPlots;

    if (this.isEnabled) {
      this.availableOptions.pop();
      this.activeOptions.pop();

      if (this.activeOptions.length > 0) {
        this.hasDropDown = this.availableOptions[this.activeOptions.length - 1].length > 7;
      }

      this.updateAllowedOptions();
      this.updateImage();
    }
  }

  @HostListener('document:click', ['$event'])
  clickout(event: { target: any; }) {
    if (!this.eRef.nativeElement.contains(event.target)) {
      this.dropdownVisible = false;
    }
  }

  updateActiveOption(level: number, newValue: string) {
    this.activeOptions[level] = newValue;
    this.updateAllowedOptions();
    this.updateImage();
  }

  updateAllowedOptions() {
    let subMap = this.data;
    let level = 0;

    for (let activeOption of this.activeOptions) {
      this.allowedOptions[level] = [];
      Object.keys(subMap).forEach(key => this.allowedOptions[level].push(key));

      if (!this.allowedOptions[level].includes(this.activeOptions[level])) {
        for (let availableOption of this.availableOptions[level]) {
          if (this.allowedOptions[level].includes(availableOption)) {
            this.activeOptions[level] = availableOption;
            break;
          }
        }
      }

      level = level + 1;
      // @ts-ignore
      subMap = subMap[activeOption];
    }
  }

  updateImage() {
    let subMap = this.data;

    for (let level of this.activeOptions) {
      // @ts-ignore
      subMap = subMap[level];
    }

    // @ts-ignore
    this.imageSource = dataPrefix + subMap["plot"];
  }

  downloadData() {
    let subMap = this.data;

    for (let level of this.activeOptions) {
      // @ts-ignore
      subMap = subMap[level];
    }

    const link = document.createElement("a");
    link.setAttribute("target", "_blank");
    // @ts-ignore
    let file: string = subMap["data"];
    link.setAttribute("href", dataPrefix + file);
    link.setAttribute("download", "test.csv");
    document.body.appendChild(link);
    link.click();
    link.remove();
  }

  collectLevels(map: {}, level: number) {
    if (this.availableOptions.length <= level) {
      this.availableOptions[level] = [];
    }

    for (let key of Object.keys(map)) {
      if (!this.availableOptions[level].includes(key)) {
        this.availableOptions[level].push(key);
      }
      if (this.activeOptions.length <= level) {
        this.activeOptions.push(key);
      }

      // @ts-ignore
      let sub = map[key];
      if (!(sub instanceof String || typeof sub == 'string')) {
        this.collectLevels(sub, level + 1);
      }
    }
  }

  moveDropdown(delta: number) {
    let dropDownLevel = this.activeOptions.length - 1;
    let currentIndex = this.allowedOptions[dropDownLevel].indexOf(this.activeOptions[dropDownLevel]);

    let newIndex = currentIndex + delta;

    if (newIndex < 0) {
      newIndex = this.allowedOptions[dropDownLevel].length - 1;
    } else if (newIndex >= this.allowedOptions[dropDownLevel].length) {
      newIndex = 0;
    }

    this.updateActiveOption(dropDownLevel, this.allowedOptions[dropDownLevel][newIndex]);
  }
}
