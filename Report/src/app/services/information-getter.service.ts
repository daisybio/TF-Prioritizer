import { Injectable } from '@angular/core';
import {INFORMATION} from "../../assets/information";

@Injectable({
  providedIn: 'root'
})
export class InformationGetterService {
  information = INFORMATION;

  constructor() { }

  getInformation() {
    return this.information;
  }
}
