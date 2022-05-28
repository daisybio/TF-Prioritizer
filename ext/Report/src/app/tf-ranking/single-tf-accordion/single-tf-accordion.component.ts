import {Component, Input, OnInit} from '@angular/core';
import {TranscriptionFactor} from "../../types/types";
import {InformationGetterService} from "../../services/information-getter.service";

@Component({
  selector: 'single-tf-accordion',
  templateUrl: './single-tf-accordion.component.html',
  styleUrls: ['./single-tf-accordion.component.css']
})
export class SingleTfAccordionComponent implements OnInit {
  @Input()
    // @ts-ignore
  transcriptionFactor: TranscriptionFactor;

  @Input()
  hasSibling: boolean = false;

  visible: boolean = false;

  information;

  constructor(private informationGetter: InformationGetterService) {
    this.information = informationGetter.getInformation();
  }

  ngOnInit(): void {
  }

  toggleVisibility() {
    this.visible = !this.visible;
  }
}
